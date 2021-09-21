"""Unit tests for pdb_structure module, including sr.PDBStructure."""

import copy
from collections import defaultdict

from Bio import SeqIO
import pytest

import schemarecomb as sr
from schemarecomb.pdb_structure import Atom, AminoAcid


@pytest.fixture
def atom_lines(bgl3_pdb_filename):
    with open(bgl3_pdb_filename) as f:
        lines = [line for line in f.read().split('\n') if line[:4] == 'ATOM'
                 and line[21] == 'A']  # chain A only.
    return lines


@pytest.fixture
def atom_list(atom_lines):
    return [Atom.from_line(line) for line in atom_lines]


@pytest.fixture
def amino_acid_list(atom_list):
    index_to_atom = defaultdict(list)
    for atom in atom_list:
        index_to_atom[atom.res_index].append(atom)
    return [AminoAcid(atoms) for atoms in index_to_atom.values()]


class TestAtom:
    @pytest.mark.parametrize(
        'atom_index,res_name,res_index,coords',
        [
            (0, 'ALA', 14, (-1.611, 17.176, 10.792)),
            (7, 'LEU', 15, (-1.711, 21.073, 14.992)),
            (412, 'HIS', 68, (10.860, 48.398, 12.372)),
        ]
    )
    def test_atom(self, atom_lines, atom_index, res_name, res_index, coords):
        atom_line = atom_lines[atom_index]
        atom = Atom.from_line(atom_line)

        assert atom.res_name == res_name
        assert atom.res_index == res_index
        assert (atom.x, atom.y, atom.z) == coords

    def test_invalid_lines(self, atom_lines):
        with pytest.raises(ValueError):
            Atom.from_line('')
        with pytest.raises(ValueError):
            Atom.from_line(atom_lines[0] + atom_lines[1])

    def test_lines(self, atom_lines):
        for line in atom_lines:
            atom = Atom.from_line(line)
            assert line == atom.to_line()


class TestAminoAcid:
    def test_init(self, amino_acid_list):
        assert amino_acid_list

    def test_invalid_init(self, atom_list):
        with pytest.raises(ValueError):
            AminoAcid([])

        # Pass in all atoms. Should have multiple different res_name values.
        with pytest.raises(ValueError):
            AminoAcid(atom_list)

        atom1 = copy.deepcopy(atom_list[0])
        atom2 = copy.deepcopy(atom_list[1])
        atom1.res_index = 12
        atom2.res_index = 13
        with pytest.raises(ValueError):
            AminoAcid([atom1, atom2])

    def test_lines(self, atom_lines):
        index_to_atom_lines = defaultdict(list)
        for line in atom_lines:
            res_index = int(line[22:26].strip()) - 1  # switch to 0 indexing
            index_to_atom_lines[res_index].append(line)
        assert index_to_atom_lines
        for input_lines in index_to_atom_lines.values():
            aa = AminoAcid.from_lines(input_lines)
            output_lines = aa.to_lines()
            assert all(i == o for i, o in zip(input_lines, output_lines))

    def test_invalid_lines(self):
        with pytest.raises(ValueError):
            AminoAcid.from_lines([])

    def test_renumbering(self, amino_acid_list):
        renum_aa_list = [copy.deepcopy(aa) for aa in amino_acid_list]
        for i, (aa, orig_aa) in enumerate(zip(renum_aa_list, amino_acid_list)):
            with pytest.raises(AttributeError):
                aa.orig_index

            # test renumbering
            aa.renumber(i)
            assert aa.index == i
            for atom in aa.atoms:
                assert atom.res_index == i
            assert aa.orig_index == orig_aa.index

            # test derenumbering
            aa.derenumber()
            true_index = orig_aa.index
            assert aa.index == true_index
            for atom in aa.atoms:
                assert atom.res_index == true_index
            with pytest.raises(AttributeError):
                aa.orig_index

    def test_json(self, amino_acid_list):
        assert amino_acid_list
        for in_aa in amino_acid_list:
            json_str = in_aa.to_json()
            assert json_str
            out_aa = AminoAcid.from_json(json_str)
            assert len(in_aa.atoms) == len(out_aa.atoms)
            assert in_aa.name == out_aa.name
            assert in_aa.index == out_aa.index

    def test_letter(self, amino_acid_list):
        for aa in amino_acid_list:
            letter = aa.letter
            assert isinstance(letter, str)
            assert len(letter) == 1

    def test_coords(self, amino_acid_list):
        for aa in amino_acid_list:
            arr = aa.coords
            assert arr.shape == (len(aa.atoms), 3)
            for atom, (x, y, z) in zip(aa.atoms, arr):
                assert atom.x == x
                assert atom.y == y
                assert atom.z == z


class TestPDBStructure:
    def invalid_init(self, bgl3_fasta_filename, amino_acid_list):
        p1 = str(list(SeqIO.parse(bgl3_fasta_filename, 'fasta'))[0].seq)
        with pytest.raises(ValueError):
            sr.PDBStructure(amino_acid_list, None, p1)

        p1_short = p1[:400]
        renum_pdb = sr.PDBStructure(amino_acid_list)
        renum_pdb.renumber(p1_short)
        amino_acids = renum_pdb.amino_acids
        unrenum = renum_pdb.unrenumbered_amino_acids
        with pytest.raises(ValueError):
            sr.PDBStructure(amino_acids, unrenum, None)

        # Empty unrenumbered_amino_acids.
        with pytest.raises(ValueError):
            sr.PDBStructure(amino_acids, [], p1_short)

        # p1 too short.
        with pytest.raises(ValueError):
            sr.PDBStructure(amino_acids, unrenum, p1_short[:200])

        # p1 too long.
        with pytest.raises(ValueError):
            sr.PDBStructure(amino_acids, unrenum, p1)

    def test_renumbering(self, bgl3_fasta_filename, amino_acid_list):
        # TODO: Probably want to do this one with multiple PDBs.
        assert amino_acid_list

        # Test sr.PDBStructure init.
        pdb = sr.PDBStructure(amino_acid_list)
        assert pdb.amino_acids
        with pytest.raises(AttributeError):
            pdb.unrenumbered_amino_acids
        with pytest.raises(AttributeError):
            pdb.renumbering_seq
        assert pdb.seq
        assert pdb.contacts

        renum_pdb = copy.deepcopy(pdb)
        p1 = str(list(SeqIO.parse(bgl3_fasta_filename, 'fasta'))[0].seq)
        renum_pdb.renumber(p1)
        assert len(pdb.amino_acids) == len(renum_pdb.amino_acids) \
            + len(renum_pdb.unrenumbered_amino_acids)
        # Technically this could be true, but not for this PDB.
        assert pdb.seq != renum_pdb.seq
        assert renum_pdb.renumbering_seq == p1
        assert pdb.contacts
        # Most important!
        assert all(aa == p1[i] for i, aa in enumerate(renum_pdb.seq)
                   if aa != '-')

        renum_pdb.derenumber()
        with pytest.raises(AttributeError):
            pdb.unrenumbered_amino_acids
        with pytest.raises(AttributeError):
            pdb.renumbering_seq
        assert pdb.seq == renum_pdb.seq
        assert len(pdb.contacts) == len(renum_pdb.contacts)
        pdb_contacts = set(pdb.contacts)
        for contact in renum_pdb.contacts:
            assert contact in pdb_contacts

        # Derenumbering when not renumbered.
        with pytest.raises(AttributeError):
            pdb.derenumber()

    def test_double_renum(self, bgl3_fasta_filename, amino_acid_list):
        """Same as test_renumbering but double renum with two reduced p1."""
        pdb = sr.PDBStructure(amino_acid_list)
        p1 = str(list(SeqIO.parse(bgl3_fasta_filename, 'fasta'))[0].seq)

        renum_pdb = copy.deepcopy(pdb)
        p1_short = p1[:400]
        renum_pdb.renumber(p1_short)
        assert renum_pdb.unrenumbered_amino_acids  # Should be occupied now.
        assert len(pdb.amino_acids) == len(renum_pdb.amino_acids) \
            + len(renum_pdb.unrenumbered_amino_acids)
        assert pdb.seq != renum_pdb.seq
        assert renum_pdb.renumbering_seq == p1_short
        assert all(aa == p1_short[i] for i, aa in enumerate(renum_pdb.seq)
                   if aa != '-')

        # Renumber again.
        p1_middle_removed = p1_short[:150] + p1_short[250:]
        renum_pdb.renumber(p1_middle_removed)
        assert renum_pdb.unrenumbered_amino_acids
        assert len(pdb.amino_acids) == len(renum_pdb.amino_acids) \
            + len(renum_pdb.unrenumbered_amino_acids)
        assert pdb.seq != renum_pdb.seq
        assert renum_pdb.renumbering_seq == p1_middle_removed
        assert all(aa == p1_middle_removed[i] for i, aa
                   in enumerate(renum_pdb.seq) if aa != '-')

        renum_pdb.derenumber()
        with pytest.raises(AttributeError):
            pdb.unrenumbered_amino_acids
        with pytest.raises(AttributeError):
            pdb.renumbering_seq
        assert pdb.seq == renum_pdb.seq
        assert len(pdb.contacts) == len(renum_pdb.contacts)
        # pdb_contacts = set(pdb.contacts)
        for contact in renum_pdb.contacts:
            assert contact in pdb.contacts

    def test_from_pdb_file(self, bgl3_pdb_filename, amino_acid_list):
        pdb_s = sr.PDBStructure.from_pdb_file(str(bgl3_pdb_filename))
        with open(bgl3_pdb_filename) as f:
            pdb_f = sr.PDBStructure.from_pdb_file(f)
        pdb_p = sr.PDBStructure.from_pdb_file(bgl3_pdb_filename)

        pdb = sr.PDBStructure(amino_acid_list)

        assert pdb.seq == pdb_s.seq
        assert pdb.seq == pdb_f.seq
        assert pdb.seq == pdb_p.seq

    def test_json(self, amino_acid_list, bgl3_fasta_filename):
        in_pdb = sr.PDBStructure(amino_acid_list)
        json_str = in_pdb.to_json()
        assert json_str
        out_pdb = sr.PDBStructure.from_json(json_str)
        assert in_pdb.seq == out_pdb.seq
        with pytest.raises(AttributeError):
            out_pdb.unrenumbered_amino_acids
        with pytest.raises(AttributeError):
            out_pdb.renumbering_seq

        p1 = str(list(SeqIO.parse(bgl3_fasta_filename, 'fasta'))[0].seq)
        p1_short = p1[:400]
        renum_pdb = sr.PDBStructure(amino_acid_list)
        renum_pdb.renumber(p1_short)
        json_str2 = renum_pdb.to_json()
        assert json_str2
        out_pdb2 = sr.PDBStructure.from_json(json_str2)
        assert renum_pdb.seq == out_pdb2.seq
        assert renum_pdb.renumbering_seq == out_pdb2.renumbering_seq
