# pdb_structure.py

"""Module for obtaining and manipulating PDB files."""

from dataclasses import dataclass
from itertools import product

from Bio import SeqIO, SeqRecord

from blast_query import blast_query
from utils import _calc_identity


@dataclass
class Atom:
    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
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
    def from_line(cls, line):
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

    def to_line(self):
        line = list('ATOM') + [' ' for _ in range(76)]
        # line[6:11] = list(str(self.serial_num).rjust(5))
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

    @property
    def coords(self):
        return (self.x, self.y, self.z)

    def d(self, atom):
        coords = zip(self.coords, atom.coords)
        return sum((c1 - c2)**2 for c1, c2 in coords) ** (1/2)


@dataclass
class AminoAcid:
    atoms: list[Atom]
    res_name: str
    res_seq_num: int

    @classmethod
    def from_lines(cls, lines):
        atoms = [Atom.from_line(line) for line in lines]
        return cls.from_atoms(atoms)

    @classmethod
    def from_atoms(cls, atoms):
        res_name = atoms[0].res_name
        assert all(a.res_name == res_name for a in atoms)
        res_seq_num = atoms[0].res_seq_num
        assert all(a.res_seq_num == res_seq_num for a in atoms)
        return cls(atoms, res_name, res_seq_num)

    def to_lines(self):
        return '\n'.join([a.to_line() for a in self.atoms])

    def d(self, res):
        return min(a1.d(a2) for a1, a2 in product(self.atoms, res.atoms))


class PDBStructure:
    """Protein crytal structure from PDB."""
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
