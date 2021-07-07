# pdb_structure.py

"""Module for obtaining and manipulating PDB files for SCHEMA-RASPP.

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
with the rest of the module, these classes index the residue number starting
from 0. As a result, PDB file reading and writing PDB must convert between
indexing. For example, the "ALA" lines in the pdb file below indicate that the
15th amino acid in the structure is alanine. When read with this module, this
alanine will be labeled with a res_seq_num of 14. This is consistent with
sequence number: if pdb_seq is the Python String representing the amino acid
sequence of the structure , this residue would be pdb_seq[14].

PDB files generally have this structure (example structure 1GNX):
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

Available classes:
    PDBStructure: Protein structure from PDB file.
    AminoAcid: Residue in PDB structure.
    Atom: Atom within a certain residue in PDB structure.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations, product
from io import BufferedIOBase, TextIOBase
from typing import Union
from urllib.request import urlopen

from Bio import pairwise2, SeqRecord
from Bio.SeqUtils import seq1
import numpy as np
from scipy.spatial import distance

from .blast_query import blast_query
from .utils import _calc_identity


@dataclass
class Atom:
    """Atom within a PDB structure.

    This class is basically a one-to-one conversion of 'ATOM' lines in PDB
    files. There are additional methods for SCHEMA-RASPP calculations.

    Input and output line formats derive from the PDB format at:
    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

    Public attributes:
        serial_num: Atom serial number.
        name: Atom name.
        alt_loc: Alternate location indicator.
        res_name: Reside name.
        chain: Chain identifier.
        res_seq_num: Residue sequence number.
        insertion_code: Code for insertions of residues.
        x: X orthogonal Angstrom coordinate.
        y: Y orthogonal Angstrom coordinate.
        z: Z orthogonal Angstrom coordinate.
        occup: Occupancy.
        temp: Temperature factor.
        segment: Segment identifier.
        element: Element symbol.
        charge: Charge.

    Public methods:
        d: Distance between this atom and another in Angstroms.
        to_line: Convert instance into PDB file ATOM line.

    Constructors:
        from_line: Build from PDB file ATOM line.
    """
    serial_num: int
    name: str
    alt_loc: str
    res_name: str
    chain: str
    res_seq_num: int
    insertion_code: int
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

        Args:
            line: ATOM line from PDB file.
        """
        serial_num = int(line[6:11].strip())
        name = line[12:16].strip()
        alt_loc = line[16]
        res_name = line[17:20].strip()
        chain = line[21]
        res_seq_num = int(line[22:26].strip()) - 1  # switch to 0 indexing
        insertion_code = line[26]
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occup = float(line[54:60].strip())
        temp = float(line[60:66].strip())
        segment = line[72:76].strip()
        element = line[76:78].strip()
        charge = line[78:80].strip()
        return cls(serial_num, name, alt_loc, res_name, chain, res_seq_num,
                   insertion_code, x, y, z, occup, temp, segment, element,
                   charge)

    def to_line(self) -> str:
        """Convert to PDB file ATOM line."""
        line = list('ATOM') + [' ' for _ in range(76)]
        line[6:11] = str(self.serial_num).rjust(5)
        if len(self.name) == 4:
            line[12:16] = self.name.ljust(4)
        else:
            line[12:16] = self.name.rjust(2) + '  '
        line[16] = self.alt_loc
        line[17:20] = self.res_name.rjust(3)
        line[21] = self.chain
        line[22:26] = str(self.res_seq_num+1).rjust(4)  # switch to 1 indexing
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

    def d(self, atom: 'Atom') -> float:
        """Distance between this atom and another in Angstroms."""
        coords1 = (self.x, self.y, self.z)
        coords2 = (atom.x, atom.y, atom.z)
        coords = zip(coords1, coords2)
        return sum((c1 - c2)**2 for c1, c2 in coords) ** (1/2)


@dataclass
class AminoAcid:
    """Amino acid within a PDB structure.

    Collection of atoms that make up an a residue in the PDB structure, plus
    additional functions for SCHEMA-RASPP calculation.

    Public attributes:
        atoms: The residue's atoms.
        res_name: Name of the residue.
        res_seq_num: Number of the residue within the structure.
        letter: res_name converted to one letter for convenience.
        coords: x, y, and z coordinates of each atom.

    Public methods:
        to_lines: Convert to section of PDB ATOM lines.
        d: Distance between closent atoms of this and another residue.
        renumber: Change the res_seq_num of the amino acid and its atoms.

    Constructors:
        from_lines: From ATOM lines in a PDB files.
        from_atoms: From list of Atoms.
    """
    atoms: list[Atom]
    res_name: str
    res_seq_num: int
    letter: str

    @property
    def coords(self) -> np.ndarray:  # shape (# atoms, 3)
        """x, y, and z coordinates in a numpy array."""
        return np.asarray([(a.x, a.y, a.z) for a in self.atoms])

    @classmethod
    def from_lines(cls, lines: Union[str, Iterable[str]]) -> 'AminoAcid':
        """Construct from ATOM lines in a PDB file."""
        if isinstance(lines, str):
            lines = lines.split('\n')
        atoms = [Atom.from_line(line) for line in lines]
        return cls.from_atoms(atoms)

    @classmethod
    def from_atoms(cls, atoms: list[Atom]) -> 'AminoAcid':
        """Construct from list of Atoms."""
        res_name = atoms[0].res_name
        assert all(a.res_name == res_name for a in atoms)
        res_seq_num = atoms[0].res_seq_num
        assert all(a.res_seq_num == res_seq_num for a in atoms)
        three_letters = res_name[0] + res_name[1:].lower()
        letter = seq1(three_letters)
        return cls(atoms, res_name, res_seq_num, letter)

    def to_lines(self) -> str:
        """Convert to section of PDB ATOM lines."""
        return '\n'.join([a.to_line() for a in self.atoms])

    def d(self, res: 'AminoAcid') -> float:
        """Distance between closest atoms of this and another residue."""
        return min(a1.d(a2) for a1, a2 in product(self.atoms, res.atoms))

    def renumber(self, new_num: int) -> None:
        """Change the res_seq_num of the amino acid and its atoms."""
        self.res_seq_num = new_num
        for a in self.atoms:
            a.res_seq_num = new_num


@dataclass
class PDBStructure:
    """Protein crytal structure from PDB.

    Container for list of amino acids.

    Public attributes:
        amino_acids: Resides in PDB structure.
        renumbered: Whether the structure has be renumbered to an alignment.
        seq: Amino acid sequence of structure.
        aligned_seq: seq aligned to a ParentAlignment.
        contacts: Indices of contacting residues. Used in SCHEMA energy.

    Public methods:
        renumber: Renumber pdb structure to match ParentAlignment.

    Constructors:
        from_parents: Initialize from aligned sequences.
        from_pdb_file: Build directly from PDB file.

    Future: This constructor's BLAST search and filter are effectively
        pointless. This constructor should probably be broken into a few
        separate constructors:
            - Input aligned PDB sequence.
            - Input PDB accession and align to parent_alignment?
            - Input parent sequences and find if any parents are PDB/closest.
            - Basically this constructor but handle PDB differences.
            - Add renumbering as a post-initialization option to enable manual
                construction of structures.
    """
    amino_acids: list[AminoAcid]
    is_renumbered: bool = False

    @property
    def seq(self) -> str:
        """Amino acid sequence of structure, only includes AAs in structure."""
        return ''.join([aa.letter for aa in self.amino_acids])

    @property
    def aligned_seq(self) -> str:
        """seq aligned to a ParentlAlignment."""
        seq_dict = {aa.res_seq_num: aa.letter for aa in self.amino_acids}
        return ''.join([seq_dict.get(i, '-') for i in range(max(seq_dict)+1)])

    @cached_property
    def contacts(self) -> list[tuple[int, int]]:
        """Lazy evaluation of contacts.

        TODO: Refactor this to be easier to use in SCHEMA-RASPP. Options:
            - Leave as is.
            - Change to pre-computed rather than lazy.
            - Change to lazy on per-AA basis.
            - Add PDBStructure method with two AA input (dict?).
            - Add contacting residues to AminoAcid object.
        """
        contacts = []
        aa_coords = [aa.coords for aa in self.amino_acids]
        combos = combinations(range(len(aa_coords)), 2)
        for ii, j in combos:
            a1 = aa_coords[ii]
            a2 = aa_coords[j]
            d = distance.cdist(a1, a2)
            if np.any(d < 4.5):
                contacts.append((ii, j))
        return contacts

    @classmethod
    def from_parents(cls,
                     parent_seqs: list[SeqRecord.SeqRecord],
                     p1_aligned: Union[str, SeqRecord.SeqRecord],
                     ) -> 'PDBStructure':
        """Construct from aligned sequences using BLAST and PDB.

        The best structure is found by using BLAST to download candidate PDB
        sequences, then the sequence with the largest minimum identity to the
        parents is selected. This sequence is aligned with the first parent to
        renumber the PDB residues to parent alignment.

        Note that the current implementation doesn't renumber PDB residues that
        aren't found in the first parent. This effectively means that the first
        parent should always be selected from the PDB database.

        Args:
            parent_seqs: Unaligned sequences from the parent alignment. Note
                that parent_seqs[0] is special: it is used to query BLAST and
                align the chosen PDB to the sequence.
            p1_aligned: parent_seqs[0] aligned to the other parents.
        """
        query_str = str(parent_seqs[0].seq)
        if isinstance(p1_aligned, SeqRecord.SeqRecord):
            p1_aligned = str(p1_aligned.seq)
        if query_str != p1_aligned.replace('-', ''):
            raise ValueError('parent_seqs[0] does not match p1_aligned.')
        pdb_srs = list(blast_query(query_str, 'pdbaa', 100))

        # Find the PDB struct with largest minimum identity to the parents.
        best_id = None  # track the best sequence
        best_min_iden = 0.0  # minimum parental identity of best sequence
        for pdb_s in pdb_srs:
            min_iden = 1.0  # minimum parental identity of current sequence
            is_minimum = True  # whether the current sequence is the best
            for par_s in parent_seqs:
                iden = _calc_identity(pdb_s, par_s)
                if iden < best_min_iden:  # sequence suboptimal, short-circuit
                    is_minimum = False
                    break
                min_iden = min(min_iden, iden)
            if is_minimum:  # never set false, must be optimal so far
                best_id = pdb_s.id
                best_min_iden = min_iden

        if best_id is None:
            raise ValueError('No best PDB found.')

        # Get the pdb_structure from rcsb.
        acc = best_id.split('|')[1]
        url = 'https://files.rcsb.org/view/' + acc + '.pdb'
        from urllib.request import Request
        request = Request(url)
        with urlopen(request) as f:
            pdb_structure = cls.from_pdb_file(f)

        pdb_structure.renumber(p1_aligned)

        return pdb_structure

    @classmethod
    def from_pdb_file(cls, f: Union[str, TextIOBase, BufferedIOBase],
                      chain: str = 'A') -> 'PDBStructure':
        """Construct from PDB file without renumbering.

        Args:
            f: file-like PDB structure.
            chain: Chain to include in constructed object.
        """
        if isinstance(f, str):
            f = open(f)
        amino_acids = []
        curr_atoms = []
        for line in f:
            if line[:4] not in (b'ATOM', 'ATOM'):
                continue
            if isinstance(line, bytes):
                line = line.decode()
            atom = Atom.from_line(line)
            if atom.chain != chain:
                continue
            if curr_atoms and curr_atoms[0].res_seq_num != atom.res_seq_num:
                aa = AminoAcid.from_atoms(curr_atoms)
                amino_acids.append(aa)
                curr_atoms = [atom]
            else:
                curr_atoms.append(atom)
        if curr_atoms:
            aa = AminoAcid.from_atoms(curr_atoms)
            amino_acids.append(aa)
        return cls(amino_acids)

    def renumber(self, p1_aligned: SeqRecord.SeqRecord) -> None:
        """Renumber pdb structure to match ParentAlignment."""
        pdb_atom_seq = ''.join([aa.letter for aa in self.amino_acids])
        p1_seq = str(p1_aligned.seq)
        aln = pairwise2.align.globalxx(p1_seq, pdb_atom_seq, gap_char='.',
                                       one_alignment_only=True)[0]
        pdb_iter = iter(self.amino_acids)
        par_index = 0  # PDB files normally index from 0.
        for i, (par_aa, pdb_aa) in enumerate(zip(aln.seqA, aln.seqB)):
            if pdb_aa != '.':
                next(pdb_iter).renumber(par_index)
            if par_aa == '.':
                # TODO: handle this more formally?
                print(f'warning: residue {i} in pdb but not par')
            else:
                par_index += 1

        self.is_renumbered = True
