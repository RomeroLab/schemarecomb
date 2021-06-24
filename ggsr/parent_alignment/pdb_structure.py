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
atoms to match SCHEMA-RASPP parent alignment. The main goal of this module is
to relieve user confusion about this process. Note that all atom/residue
numbering in this module is indexed starting from 1 in order to match the PDB
format. Since Python is indexed starting from 0, you must be cognizant of how
this affects indexing vs labeling. For example, the N in ALA below will have
serial number 1, but will be ala15.atoms[0] if "ala15" is the AminoAcid
instance representing the ALA at position 15 in the structure. Similarly, the
"15" in "ala15" states that this residue is the 15th element in the structure.
However, if pdb_seq was the Python String representing the amino acid sequence
of the structure below, this residue would be pdb_seq[14].

PDB files generally have this structure (example structure is 1GNX):
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

from collections.abs import Iterable
from dataclasses import dataclass
from itertools import product
from typing import Union

from Bio import SeqIO, SeqRecord

from blast_query import blast_query
from utils import _calc_identity


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
        res_seq_num = int(line[22:26].strip())
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
        line[22:26] = str(self.res_seq_num).rjust(4)
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

    Public methods:
        to_lines: Convert to section of PDB ATOM lines.
        d: Distance between closent atoms of this and another residue.

    Constructors:
        from_lines: From ATOM lines in a PDB files.
        from_atoms: From list of Atoms.
    """
    atoms: list[Atom]
    res_name: str
    res_seq_num: int

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
        return cls(atoms, res_name, res_seq_num)

    def to_lines(self) -> str:
        """Convert to section of PDB ATOM lines."""
        return '\n'.join([a.to_line() for a in self.atoms])

    def d(self, res: 'AminoAcid') -> float:
        """Distance between closest atoms of this and another residue."""
        return min(a1.d(a2) for a1, a2 in product(self.atoms, res.atoms))


class PDBStructure:
    """Protein crytal structure from PDB."""
    # TODO: File this class out further.
    def __init__(self):
        pass

    @classmethod
    def from_sr_list(cls, fn: str, in_sr):
        pdb_srs = list(SeqIO.parse(fn, 'fasta'))
        chosen_sr = max(pdb_srs, key=lambda x: _calc_identity(in_sr, x))
        print(chosen_sr)

    @classmethod
    def from_query_sr(cls, in_sr: SeqRecord.SeqRecord):
        seq_str = str(in_sr.seq)
        pdb_srs = list(blast_query(seq_str, 'pdbaa'))
        fn = 'temp.fasta'
        SeqIO.write(pdb_srs, fn, 'fasta')
        cls.from_sr_list(fn, in_sr)

    @classmethod
    def from_pdb_file(cls, fn):
        amino_acids = []
        with open(fn) as f:
            curr_atoms = []
            for line in f:
                if line[:4] != 'ATOM':
                    continue
                atom = Atom.from_line(line)
                if curr_atoms \
                   and curr_atoms[0].res_seq_num != atom.res_seq_num:
                    aa = AminoAcid.from_atoms(curr_atoms)
                    amino_acids.append(aa)
                    curr_atoms = [atom]
                else:
                    curr_atoms.append(atom)
        return amino_acids


if __name__ == '__main__':
    # sr = list(SeqIO.parse('../../tests/bgl3_sample/bgl3_sequences.fasta',
    #                       'fasta'))[0]
    # pdb = PDBStructure.from_query_sr(sr)
    # pdb = PDBStructure.from_sr_list('temp.fasta', sr)
    amino_acids = PDBStructure.from_pdb_file('1gnx.pdb')
    print(len(amino_acids))
    print(amino_acids[0].to_lines())
