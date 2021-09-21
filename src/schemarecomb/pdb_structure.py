# pdb_structure.py

"""Parsing and manipulation of Protein Data Bank structure files.

This module provides the definition of :class:`schemarecomb.PDBStructure` and
the accessory classes :class:`AminoAcid` and :class:`Atom` that represent the
eponymous entities within a PDB structure.

PDBStructures contain a list of AminoAcid objects, which act as containers for
Atom objects that are read from "ATOM" lines in a PDB structure file.

Parsing and modifying protein structure data is necessary for recombinant
library design because SCHEMA energy calculations require the distance between
amino acids within the protein.

These structures are obtained from the Protein Data Bank at
https://www.rcsb.org and loaded as Python objects. Nearly all information is
discarded except for ATOM lines, which specify data about atoms within the
structure.

The most confusing part about PDB structure manipulation is the renumbering of
atoms to match SCHEMA-RASPP parent alignment. Relieving user confusion about
this process is a primary goal of this module. Note that atom/residue indexing
in PDB files begins at 1, while Python's indexing starts at 0. For consistency
with the rest of the package, these classes index the residue number starting
from 0. As a result, PDB file reading and writing PDB must convert between
indexing. For example, the "ALA" lines in the pdb file below indicate that the
15th amino acid in the structure is alanine. When read with this module, this
alanine will be labeled with an index of 14. This is consistent with
sequence number: if pdb_seq is the Python String representing the amino acid
sequence of the structure , this residue would be pdb_seq[14].

PDB files generally have this structure (example structure 1GNX)::

 ...<other structure data>...
 ATOM      1  N   ALA A  15      -1.611  17.176  10.792  1.00 36.46           N
 ATOM      2  CA  ALA A  15      -1.871  18.610  11.107  1.00 36.85           C
 ATOM      3  C   ALA A  15      -2.021  18.795  12.611  1.00 36.41           C
 ATOM      4  O   ALA A  15      -2.983  18.321  13.215  1.00 38.36           O
 ATOM      5  CB  ALA A  15      -3.131  19.081  10.392  1.00 35.10           C
 ATOM      6  N   LEU A  16      -1.064  19.496  13.206  1.00 34.22           N
 ATOM      7  CA  LEU A  16      -1.061  19.738  14.642  1.00 29.97           C
 ATOM      8  C   LEU A  16      -1.711  21.073  14.992  1.00 30.29           C
 ATOM      9  O   LEU A  16      -1.462  22.089  14.341  1.00 30.34           O
 ATOM     10  CB  LEU A  16       0.380  19.716  15.152  1.00 26.33           C
 ATOM     11  CG  LEU A  16       1.228  18.548  14.639  1.00 24.12           C
 ATOM     12  CD1 LEU A  16       2.681  18.761  15.026  1.00 22.75           C
 ATOM     13  CD2 LEU A  16       0.698  17.230  15.195  1.00 23.37           C
 ATOM     14  N   THR A  17      -2.541  21.066  16.028  1.00 30.68           N
 ATOM     15  CA  THR A  17      -3.217  22.278  16.472  1.00 29.15           C
 ATOM     16  C   THR A  17      -2.576  22.784  17.756  1.00 28.03           C
 ATOM     17  O   THR A  17      -2.337  22.012  18.683  1.00 28.04           O
 ATOM     18  CB  THR A  17      -4.716  22.014  16.733  1.00 30.86           C
 ATOM     19  OG1 THR A  17      -5.357  21.666  15.501  1.00 32.64           O
 ATOM     20  CG2 THR A  17      -5.385  23.246  17.319  1.00 30.50           C
 ...<other atoms>...
 ...<other structure data>...

"""

from dataclasses import dataclass
from functools import cached_property
from itertools import combinations
import json
from pathlib import Path
from typing import Optional, TextIO, Union

from Bio import pairwise2
from Bio.SeqUtils import seq1
import numpy as np
from scipy.spatial import distance


@dataclass
class Atom:
    """Atom within a PDB structure.

    This class is basically a one-to-one conversion of 'ATOM' lines in `PDB
    files <https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro
    .html>`_. The primary constructor is Atom.from_line() to allow for easy
    parsing of PDB files. For schemarecomb, the important attributes are
    res_name, res_index, x, y, and z.

    Initalization parameters and attributes are equivalent for this class.

    Parameters:
        serial_num: Atom serial number.
        name: Atom name.
        alt_loc: Alternate location indicator.
        res_name: Residue name.
        chain: Chain identifier.
        res_index: Residue index within the structure..
        insertion_code: Code for insertions of residues.
        x: X orthogonal Angstrom coordinate.
        y: Y orthogonal Angstrom coordinate.
        z: Z orthogonal Angstrom coordinate.
        occup: Occupancy.
        temp: Temperature factor.
        segment: Segment identifier.
        element: Element symbol.
        charge: Charge.

    """
    serial_num: int
    name: str
    alt_loc: str
    res_name: str
    chain: str
    res_index: int
    insertion_code: str
    x: float
    y: float
    z: float
    occup: float
    temp: float
    segment: str
    element: str
    charge: str

    @classmethod
    def from_line(cls, line: str) -> 'Atom':
        """Construct Atom from PDB ATOM line.

        Parameters:
            line: ATOM line from PDB file.

        Returns:
            Atom constructed from line.

        Raises:
            ValueError: If line has more than one newline, a newline not at the
                end of the line, first characters that are not "ATOM", or more
                than one occurence of ATOM, or if an Atom cannot otherwise
                be constructed from the line.

        """
        num_newlines = line.count('\n')
        if num_newlines > 1 or (num_newlines == 1 and line[-1] != '\n'):
            raise ValueError('More than one line passed into constructor.')

        if line[:4] != 'ATOM' or line.count('ATOM') != 1:
            raise ValueError('Invalid ATOM line passed into ATOM.from_line():'
                             f' {line}')

        try:
            serial_num = int(line[6:11].strip())
            name = line[12:16].strip()
            alt_loc = line[16]
            res_name = line[17:20].strip()
            chain = line[21]
            res_index = int(line[22:26].strip()) - 1  # switch to 0 indexing
            insertion_code = line[26]
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            occup = float(line[54:60].strip())
            temp = float(line[60:66].strip())
            segment = line[72:76].strip()
            element = line[76:78].strip()
            charge = line[78:80].strip()
        except ValueError:
            raise ValueError('Invalid ATOM line passed into ATOM.from_line():'
                             f' {line}')

        return cls(serial_num, name, alt_loc, res_name, chain, res_index,
                   insertion_code, x, y, z, occup, temp, segment, element,
                   charge)

    def to_line(self) -> str:
        """Convert to PDB file ATOM line.

        Returns:
            str of ATOM line in PDB file format.

        """
        line = list('ATOM') + [' ' for _ in range(75)]
        line[6:11] = str(self.serial_num).rjust(5)
        # Atom names are complicated, see the "misalignment" section of the
        # PDB format link, found in the Atom class docstring. We take a
        # technically incorrect, but more common approach.
        line[13:16] = self.name.ljust(4)
        line[16] = self.alt_loc
        line[17:20] = self.res_name.rjust(3)
        line[21] = self.chain
        line[22:26] = str(self.res_index+1).rjust(4)  # switch to 1 indexing
        line[26] = self.insertion_code
        line[30:38] = f'{self.x:.3f}'.rjust(8)
        line[38:46] = f'{self.y:.3f}'.rjust(8)
        line[46:54] = f'{self.z:.3f}'.rjust(8)
        line[54:60] = f'{self.occup:.2f}'.rjust(6)
        line[60:66] = f'{self.temp:.2f}'.rjust(6)
        line[72:76] = self.segment.ljust(4)
        line[76:78] = self.element.rjust(2)
        line[78:80] = self.charge.rjust(2)
        return ''.join(line)


class AminoAcid:
    """Amino acid within a PDB structure.

    This class is a collection of :class:`Atom` objects that make up a residue
    in a :class:`schemarecomb.PDBStructure`, plus additional functions for
    schemarecomb calculation.

    Parameters:
        atoms: Atoms within the AminoAcid. All Atom objects in list must have
            the same res_name and res_index. res_index must be indexed from 0,
            not 1, but this should have happened on Atom initialization.

    Attributes:
        atoms (list[Atom]): The residue's atoms.
        name (str): Name of the residue, e.g. ALA, LEU.
        index (int): Index of the residue within the structure.
        orig_index (int): Original struct_index. Not set unless the AminoAcid
            is renumbered.
        letter (str): name converted to one letter code.
        coords (np.ndarray): x, y, and z coordinates of each atom in a Nx3
            numpy ndarray, where N is the number of atoms.

    Raises:
        ValueError: During initialization, if the list is empty or atoms do not
            all have the same res_name and res_index.
        AttributeError: If AminoAcid has not been renumbered and orig_index is
            accessed.

    """
    def __init__(self, atoms: list[Atom]):
        if not atoms:
            raise ValueError('atoms must not be an empty list.')

        self.atoms = atoms

        res_name = atoms[0].res_name
        if any(a.res_name != res_name for a in atoms):
            raise ValueError('Atoms do not all have the same res_name.')
        self.name = res_name

        res_index = atoms[0].res_index
        if any(a.res_index != res_index for a in atoms):
            raise ValueError('Atoms do not all have the same res_index.')
        self.index = res_index

    @property
    def letter(self) -> str:
        three_letters = self.name[0] + self.name[1:].lower()
        return seq1(three_letters)

    @cached_property
    def coords(self) -> np.ndarray:  # shape (# atoms, 3)
        return np.asarray([(a.x, a.y, a.z) for a in self.atoms])

    def renumber(self, new_index: int) -> None:
        """Change the structure index of the amino acid and its atoms."""
        self.orig_index = self.index
        self.index = new_index
        for a in self.atoms:
            a.res_index = new_index

    def derenumber(self) -> None:
        """Revert changes made by the renumber method."""
        old_index = self.orig_index
        del self.orig_index
        self.index = old_index
        for a in self.atoms:
            a.res_index = old_index

    @classmethod
    def from_lines(cls, lines: list[str]) -> 'AminoAcid':
        """Construct AminoAcid from ATOM lines in a PDB file.

        Parameters:
            lines: ATOM lines in a PDB file.

        Returns:
            AminoAcid constructed from lines.

        """
        atoms = [Atom.from_line(line) for line in lines]
        return cls(atoms)

    def to_lines(self) -> list[str]:
        """Convert to section of PDB ATOM lines.

        Returns:
            PDB format line for each of the AminoAcid's atoms.

        Note:
            This method does not save the orig_index attribute.

        Example:
         >>> from schemarecomb.pdb_structure import AminoAcid
         >>> pdb_fn = 'tests/fixtures/bgl3_full/1GNX.pdb'
         >>>
         >>> # Get all "ATOM" lines with index 15 into in_lines.
         >>> in_lines = []
         >>> with open(pdb_fn) as f:
         ...     for line in f:
         ...         if line[:4] != 'ATOM':
         ...             continue
         ...         res_index = int(line[22:26].strip())
         ...         # 1GNX starts at 15, see pdb_structure module docs.
         ...         if res_index == 15:
         ...             in_lines.append(line)
         ...         else:
         ...             break
         ...
         >>> # Make AminoAcid out of in_lines. It has index 14 because we index
         >>> # starting from 0, instead of 1 in the PDB file.
         >>> aa = AminoAcid.from_lines(in_lines)
         >>> aa.index == 14
         True
         >>>
         >>> # to_lines method should produce same as in_lines.
         >>> out_lines = aa.to_lines()
         >>> all(i.strip() == o.strip() for i, o in zip(in_lines, out_lines))
         True

        """
        return [a.to_line() for a in self.atoms]

    @classmethod
    def from_json(cls, in_json: str) -> 'AminoAcid':
        """JSON from_lines wrapper that includes orig_index.

        See the from_json method for an example.

        Parameters:
            in_json: JSON-formatted input string of a list that contains a list
                of ATOM lines and orig_index, the latter of which may be null.
                Usually generated by the to_json method.

        Returns:
            AminoAcid constucted from input JSON string.

        """
        lines, orig_index = json.loads(in_json)
        obj = cls.from_lines(lines)
        if orig_index is not None:
            obj.orig_index = orig_index
        return obj

    def to_json(self) -> str:
        """JSON to_lines wrapper that includes orig_index.

        Returns:
            JSON string representing the AminoAcid.

        Example:
         >>> from schemarecomb import PDBStructure
         >>> from schemarecomb.pdb_structure import AminoAcid
         >>>
         >>> # Get the first amino acid in the PDB file.
         >>> pdb_fn = 'tests/fixtures/bgl3_full/1GNX.pdb'
         >>> pdb = PDBStructure.from_pdb_file(pdb_fn)
         >>> aa1 = pdb.amino_acids[0]
         >>>
         >>> # Make a JSON string representing aa1.
         >>> json_str = aa1.to_json()
         >>>
         >>> # Make a new AminoAcid from the JSON.
         >>> aa2 = AminoAcid.from_json(json_str)
         >>>
         >>> # aa1 and aa2 should be the same.
         >>> aa1.name == aa2.name
         True
         >>> aa1.index == aa2.index
         True

        """
        orig_index: Optional[int]
        try:
            orig_index = self.orig_index
        except AttributeError:
            orig_index = None
        return json.dumps([self.to_lines(), orig_index])


class _PDBStructure:
    """Structure from the Protein Data Bank. Used in downstream calculations.

    This class is a Python representation of a PDB file, e.g. `1GNX.pdb
    <https://files.rcsb.org/download/1GNX.pdb>`_, which hold protein structure
    data and are downloaded from https://www.rcsb.org. This download may be
    done automatically with the :meth:`schemarecomb.ParentSequences.get_PDB`
    method or manually. In the latter case, PDBStructure object must also be
    created manually, e.g. using the :meth:`from_pdb_file` class method.

    See :mod:`~schemarecomb.pdb_structure` for lower-level information about
    this class, including the PDB file format and the helper classes
    :class:`~schemarecomb.pdb_structure.AminoAcid` and
    :class:`~schemarecomb.pdb_structure.Atom`.

    Parameters:
        amino_acids: Amino acids in the PDB structure.
        unrenumbered_amino_acids: The amino acids not included in renumbering.
            None if PDBStructure is not renumbered. Must be None if
            renumbering_seq is None.
        renumbering_seq: Sequence used to renumber structure. None if
            PDBStructure is not renumbered. Must not None if
            unrenumbered_amino_acids is None.

    Attributes:
        amino_acids (list[AminoAcid]): Residues in the PDB structure. There is
            no guarantee that the AminoAcid objects are ordered by index.
        unrenumbered_amino_acids (list[AminoAcid]): Original amino_acids not
            included in the renumbering. Present if and only if structure is
            renumbered.
        renumbering_seq (str): Sequence that was used to renumber structure.
            Present if and only if structure is renumbered.
        seq (str): Amino acid sequence of structure, including gaps based on
            the indices of the amino_acids. Note that this attribute will
            change if renumbering occurs.
        contacts (list[tuple[int, int]]): Indices of contacting residues, where
            a pair of residues are contacting if the shortest distance between
            them is less than 4.5 angstroms.

    Raises:
        ValueError: During initialization, if input amino_acids have non-unique
            indicies, unrenumbered_amino_acids is an empty list,
            renumbering_seq does not match amino_acids, or one but not both of
            unrenumbered_amino_acids and renumbering_seq is None.
        AttributeError: If PDBStructure is not renumbered and
            unrenumbered_amino_acids or renumbering_seq is accessed.

    """

    def __init__(
        self,
        amino_acids: list[AminoAcid],
        unrenumbered_amino_acids: Optional[list[AminoAcid]] = None,
        renumbering_seq: Optional[str] = None
    ):
        # Ensure that all amino_acid have unique indices.
        indices = set()
        for aa in amino_acids:
            index = aa.index
            if index in indices:
                raise ValueError('AminoAcid objects in a PDBStructure must '
                                 'have unique indices.')
            indices.add(index)

        self.amino_acids = amino_acids
        self._contacts: list[tuple[int, int]]

        # If PDBStructure is already aligned, both unrenumberd_amino_acids and
        # renumbering_seq must both be not None.
        if unrenumbered_amino_acids is not None \
           and renumbering_seq is not None:
            if any(aa.letter != renumbering_seq[aa.index] for aa
                   in amino_acids):
                raise ValueError(f'renumbering seq "{renumbering_seq}" is '
                                 'invalid with input amino_acids.')

            self.unrenumbered_amino_acids = unrenumbered_amino_acids
            self.renumbering_seq = renumbering_seq
        elif unrenumbered_amino_acids is not None \
                or renumbering_seq is not None:
            raise ValueError('Unrenumbered_amino_acids and renumbering_seq '
                             'must both be None if one is None.')

    @property
    def seq(self) -> str:
        seq_dict = {aa.index: aa.letter for aa in self.amino_acids}
        return ''.join([seq_dict.get(i, '-') for i in range(max(seq_dict)+1)])

    @property
    def contacts(self) -> list[tuple[int, int]]:
        try:
            # Cached contacts.
            return self._contacts
        except AttributeError:
            contacts = []
            aa_coords = {aa.index: aa.coords for aa in self.amino_acids}
            combos = combinations(aa_coords, 2)
            for ii, j in combos:
                a1 = aa_coords[ii]
                a2 = aa_coords[j]
                d = distance.cdist(a1, a2)
                if np.any(d < 4.5):
                    contacts.append((ii, j))
                    contacts.append((j, ii))
            self._contacts = contacts
            return contacts

    def renumber(self, p0_aligned: str) -> None:
        """Renumber pdb structure to match a ParentSequences.

        Parameters:
            p0_aligned: Sequence to align to. Usually the first parent from
                a ParentSequences.

        """
        # If already renumbered, derenumber to get all amino_acids back.
        if hasattr(self, 'unrenumbered_amino_acids'):
            self.derenumber()

        # Will put each AminoAcid into one of these.
        renumbered_amino_acids = []
        unrenumbered_amino_acids = []

        # Order self.amino_acids.
        index_to_amino_acids = {aa.index: aa for aa in self.amino_acids}
        ordered_amino_acids = [index_to_amino_acids[i] for i
                               in sorted(index_to_amino_acids)]

        # Align p0_aligned and the PDBStructure, using '.' for gaps.
        pdb_seq = ''.join([aa.letter for aa in ordered_amino_acids])
        aln = pairwise2.align.globalxx(p0_aligned, pdb_seq, gap_char='.',
                                       one_alignment_only=True)[0]
        aligned_p1_seq = aln.seqA
        aligned_pdb_seq = aln.seqB

        # Simultaneously iterate through ordered amino_acids and the
        # alignment. Every time we see a non-gap character in aligned_pdb_seq,
        # the next element in pdb_iter is the AminoAcid corresponding to that
        # character. If p1 character at the same position in the alignment is
        # also not a gap, we renumber the AminoAcid to the parent_index and
        # include it.
        pdb_iter = iter(ordered_amino_acids)
        parent_index = 0  # Tracks the index of the p1.

        # There are three cases for each iteration of alignment_iter,
        # corresponding to whether each of par_aa and pdb_aa are '.' or not,
        # where '.' denotes a gap in that sequence. The fourth case is both
        # par_aa and pdb_aa are '.', which should not happen with the pairwise
        # alignment settings.
        alignment_iter = zip(aligned_p1_seq, aligned_pdb_seq)
        for par_aa, pdb_aa in alignment_iter:
            if pdb_aa != '.':
                candidate_aa = next(pdb_iter)
                if par_aa != '.':
                    # Aligned residue. Renumber to parent and include.
                    candidate_aa.renumber(parent_index)
                    renumbered_amino_acids.append(candidate_aa)
                    parent_index += 1
                else:
                    # In pdb but not p1, include in unrenumbered_amino_acids.
                    # print(f'warning: residue {i} in pdb but not parent.')
                    unrenumbered_amino_acids.append(candidate_aa)
            elif par_aa != '.':
                # Not in pdb, but in p1, increment parent_index.
                parent_index += 1

        self.amino_acids = renumbered_amino_acids
        self.unrenumbered_amino_acids = unrenumbered_amino_acids
        self.renumbering_seq = p0_aligned
        if hasattr(self, '_contacts'):
            del self._contacts

    def derenumber(self) -> None:
        """Revert any changes made by the renumber method back to original.

        Raises:
            AttributeError: If PDBStructure is not renumbered.

        """
        try:
            unrenumbered_amino_acids = self.unrenumbered_amino_acids
            del self.unrenumbered_amino_acids
        except AttributeError:
            raise AttributeError('PDBStructure has not been renumbered.')

        for aa in self.amino_acids:
            aa.derenumber()
        amino_acids = self.amino_acids + unrenumbered_amino_acids
        self.amino_acids = amino_acids
        del self.renumbering_seq
        if hasattr(self, '_contacts'):
            del self._contacts

    @classmethod
    def from_pdb_file(
        cls,
        f: Union[str, Path, TextIO],
        chain: str = 'A'
    ) -> '_PDBStructure':
        """Construct from PDB file without renumbering.

        Reference the :mod:`~schemarecomb.pdb_structure` module to see the
        structure of a PDB file.

        Parameters:
            f: Filename of PDB structure file or file-like PDB structure.
            chain: Chain to include in constructed object.

        Returns:
            PDBStructure initialized from PDB file.

        """
        if isinstance(f, str) or isinstance(f, Path):
            f = open(f)
            close_later = True
        else:
            close_later = False

        amino_acids = []
        curr_atoms: list[Atom] = []
        for line in f:
            if line[:4] not in (b'ATOM', 'ATOM'):
                continue
            if isinstance(line, bytes):
                line = line.decode()
            atom = Atom.from_line(line)
            if atom.chain != chain:
                continue
            if curr_atoms and curr_atoms[0].res_index != atom.res_index:
                aa = AminoAcid(curr_atoms)
                amino_acids.append(aa)
                curr_atoms = [atom]
            else:
                curr_atoms.append(atom)

        if close_later:
            f.close()

        if curr_atoms:
            aa = AminoAcid(curr_atoms)
            amino_acids.append(aa)
        return cls(amino_acids)

    def to_json(self) -> str:
        """Convert structure to JSON string.

        See the from_json method for an example.

        Returns:
            JSON string representing the PDBStructure.

        """
        aa_jsons = [aa.to_json() for aa in self.amino_acids]
        unre_jsons: Optional[list[str]]
        try:
            unre_jsons = [aa.to_json() for aa in self.unrenumbered_amino_acids]
        except AttributeError:
            unre_jsons = None
        try:
            renum_seq: Optional[str] = self.renumbering_seq
        except AttributeError:
            renum_seq = None
        return json.dumps([aa_jsons, unre_jsons, renum_seq])

    @classmethod
    def from_json(cls, in_json: str) -> '_PDBStructure':
        """Construct from JSON string.

        Parameters:
            in_json: JSON string representing a PDBStructure instance. Usually
                generated by the to_json method.

        Returns:
            PDBStructure constructed from JSON string.

        Example::

         >>> from schemarecomb import PDBStructure
         >>> from schemarecomb.pdb_structure import AminoAcid
         >>>
         >>> # Get a PDBStructure from a PDB file.
         >>> pdb_fn = 'tests/fixtures/bgl3_full/1GNX.pdb'
         >>> pdb1 = PDBStructure.from_pdb_file(pdb_fn)
         >>>
         >>> # Makes a temporary directory that can be cleaned up later. You
         >>> # can ignore this and just use a string for your filename.
         >>> tmpdir = getfixture('tmpdir')
         >>> json_filename = tmpdir / 'pdb_structure.json'
         >>>
         >>> # Convert pdb into a JSON string and save it.
         >>> json_str = pdb1.to_json()
         >>> with open(json_filename, 'w') as f:
         ...     f.write(json_str)
         ...
         305616
         >>> # Load the JSON string and make a new PDBStructure.
         >>> with open(json_filename) as f:
         ...     in_str = f.read()
         ...
         >>> pdb2 = PDBStructure.from_json(in_str)
         >>>
         >>> # pdb1 and pdb2 are the same.
         >>> aas1 = pdb1.amino_acids
         >>> aas2 = pdb2.amino_acids
         >>> aa_zip = list(zip(aas1, aas2))
         >>> all(aa1.name == aa2.name for aa1, aa2 in aa_zip)
         True
         >>> all(aa1.index == aa2.index for aa1, aa2 in aa_zip)
         True
         >>> all(len(aa1.atoms) == len(aa2.atoms) for aa1, aa2 in aa_zip)
         True

        """
        aa_jsons, unre_jsons, renumbering_seq = json.loads(in_json)
        amino_acids = [AminoAcid.from_json(aa_json) for aa_json in aa_jsons]
        unrenumbered: Optional[list[AminoAcid]]
        if unre_jsons is None:
            unrenumbered = None
        else:
            unrenumbered = [AminoAcid.from_json(aa_json) for aa_json
                            in unre_jsons]
        return cls(amino_acids, unrenumbered, renumbering_seq)
