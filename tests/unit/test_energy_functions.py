import pytest

# from ggrecomb import ParentSequences
# from ggrecomb import PDBStructure
import ggrecomb as sr
from ggrecomb.energy_functions import SCHEMA


@pytest.fixture
def bgl3_trunc_dir(fixture_dir):
    return fixture_dir / 'bgl3_50-aa'


@pytest.fixture
def bgl3_trunc_parents(bgl3_trunc_dir):
    pdb = sr.PDBStructure.from_pdb_file(bgl3_trunc_dir / '1GNX_trunc.pdb')
    parents_fn = bgl3_trunc_dir / 'bgl3_trunc_aln.fasta'
    parents = sr.ParentSequences.from_fasta(parents_fn, pdb_structure=pdb,
                                            prealigned=True)
    return parents


def test_init(bgl3_trunc_parents, bgl3_records, bgl3_records_aln):
    # Good: aligned and has pdb_structure
    assert hasattr(bgl3_trunc_parents, 'alignment')
    assert hasattr(bgl3_trunc_parents, 'pdb_structure')
    schema = SCHEMA(bgl3_trunc_parents)
    assert hasattr(schema, 'E_matrix')

    # No alignment.
    parents = sr.ParentSequences(bgl3_records)
    assert not hasattr(parents, 'alignment')
    with pytest.raises(ValueError):
        SCHEMA(parents)

    # Alignment but no pdb_structure.
    aln_parents = sr.ParentSequences(bgl3_records_aln, prealigned=True)
    assert hasattr(aln_parents, 'alignment')
    assert not hasattr(aln_parents, 'pdb_structure')
    with pytest.raises(ValueError):
        SCHEMA(aln_parents)


def eq_6(start, end, parents):
    num_chims = len(parents.records) ** 2
    aln = parents.alignment
    sequences = list(zip(*aln))
    contacts = set(parents.pdb_structure.contacts)

    def chimera_e(seqr, seqt, blkr_start, blkt_start):
        # Inside (r and t sums) of eq 6.
        if seqr is seqt:
            return 0.0
        energy = 0.0
        for r in range(blkr_start, blkt_start):
            for t in range(blkt_start, len(aln)):
                if (r, t) not in contacts:
                    continue
                aasr = aln[r]
                aast = aln[t]
                if (seqr[r], seqt[t]) not in zip(aasr, aast):
                    energy += 1.0
        return energy

    manual_energy_sum = 0.0
    for seqr in sequences:
        for seqt in sequences:
            manual_energy_sum += chimera_e(seqr, seqt, start, end)
    # Don't need to do subtraction part of eq. 6 b/c it will always be 0 in
    # SCHEMA.
    return manual_energy_sum / num_chims


def test_block_energy(bgl3_trunc_parents):
    # Comparing SCHEMA object to full Endelman et al 2004 eq 6 calculation.
    schema = SCHEMA(bgl3_trunc_parents)

    energies = {}

    for start in range(len(bgl3_trunc_parents.alignment) - 1):
        for end in range(start+1, len(bgl3_trunc_parents.alignment)):
            man_energy = eq_6(start, end, bgl3_trunc_parents)
            assert abs(schema.block_energy(start, end) - man_energy) < 10**-6
            energies[(start, end)] = man_energy

    for start in range(len(bgl3_trunc_parents.alignment) - 1):
        curr_e = energies[(start, start+1)]
        for end in range(start+2, len(bgl3_trunc_parents.alignment)):
            curr_e += schema.increment_block_energy(start, end)
            assert abs(curr_e - energies[(start, end)]) < 10**-6
