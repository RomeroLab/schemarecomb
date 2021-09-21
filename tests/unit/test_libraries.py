from collections import deque
from decimal import Decimal
from itertools import product
from operator import attrgetter
import random

import pytest

import schemarecomb as sr


@pytest.fixture
def schema(bgl3_parent_alignment):
    return sr.energy_functions.SCHEMA(bgl3_parent_alignment)


@pytest.fixture
def bsaI():
    return sr.restriction_enzymes.RestrictionEnzyme.from_name('BsaI-HFv2')


@pytest.fixture
def breakpoints(bgl3_parent_alignment, AA_C31):
    return sr.breakpoints.calculate_breakpoints(bgl3_parent_alignment, AA_C31)


@pytest.fixture
def mr_cache(bgl3_parent_alignment, breakpoints):
    return sr.libraries.MutationRateCache.from_parents(bgl3_parent_alignment,
                                                       breakpoints)


@pytest.fixture
def lib_config(schema, bsaI, mr_cache):
    return sr.libraries.LibraryConfig(schema, bsaI, 0.99, mr_cache)


def generate_bp_sets(breakpoints, max_ind, num_bps, min_blk_len, max_blk_len):
    bp_sets = []
    sorted_bps = sorted(breakpoints, key=attrgetter('position'))

    stack = deque((bp,) for bp in sorted_bps
                  if min_blk_len < bp.position < max_blk_len)
    while stack:
        bps = stack.pop()

        if len(bps) == num_bps:
            if min_blk_len < max_ind - bps[-1].position < max_blk_len:
                bp_sets.append(bps)
            continue

        indi = bps[-1].position + min_blk_len
        indf = bps[-1].position + max_blk_len
        for bp in sorted_bps:
            if indi < bp.position < indf:
                stack.append(bps + (bp,))

    return bp_sets


def test_average_mutations(breakpoints, bgl3_parent_alignment):
    # compare to brute force
    max_ind = len(bgl3_parent_alignment.alignment)
    bp_sets = generate_bp_sets(breakpoints, max_ind, 3, 110, 150)

    aln_seqs = [''.join(seq) for seq in zip(*bgl3_parent_alignment.alignment)]

    bp_sample = random.sample(bp_sets, 5)

    for bp_set in bp_sample:
        m_shortcut = sr.libraries.average_mutations(bp_set,
                                                    bgl3_parent_alignment)

        # brute force
        m_tot = 0.0
        block_inds = sr.breakpoints.block_indices(bp_set, max_ind)
        for chim_parents in product(aln_seqs, repeat=len(block_inds)):
            chim_list = [par[start:end] for par, (start, end)
                         in zip(chim_parents, block_inds)]
            chim_seq = ''.join(chim_list)

            mut_nums = [sr.libraries._sequence_mutations(chim_seq, par)
                        for par in aln_seqs]

            m_tot += min(mut_nums)

        m_brute = Decimal(m_tot) / Decimal(len(aln_seqs) ** len(block_inds))

        assert abs(m_shortcut - m_brute) < 10**-6


def test_mr_cache(breakpoints, bgl3_parent_alignment, mr_cache):
    # all calculated vs cache
    max_ind = len(bgl3_parent_alignment.alignment)
    bp_sets = generate_bp_sets(breakpoints, max_ind, 3, 110, 150)

    for i, bp_set in enumerate(bp_sets):
        m_no_cache = sr.libraries.average_mutations(bp_set,
                                                    bgl3_parent_alignment)

        m_cache = sr.libraries.average_mutations(bp_set, bgl3_parent_alignment,
                                                 mr_cache)

        assert abs(m_cache - m_no_cache) < 10**6

    # make sure there were cache hits
    # Possibly sampling error, but more likely a code error.
    assert len(bp_sets) != len(mr_cache)


def test_config(lib_config, breakpoints):
    parents = lib_config.energy_function.parents
    max_ind = len(parents.alignment)
    bp_sets = generate_bp_sets(breakpoints, max_ind, 4, 90, 110)

    bp_sample = random.sample(bp_sets, 5)

    for bp_set in bp_sample:
        assert sr.Library.calc_from_config(bp_set, Decimal(1.0), lib_config)

    lib_config_no_cache = sr.libraries.LibraryConfig(
        lib_config.energy_function,
        lib_config.gg_enzyme,
        lib_config.gg_threshold,
        None
    )

    for bp_set in bp_sample:
        assert sr.Library.calc_from_config(bp_set, Decimal(1.0),
                                           lib_config_no_cache)


def test_json(lib_config, breakpoints):
    parents = lib_config.energy_function.parents
    max_ind = len(parents.alignment)
    bp_sets = generate_bp_sets(breakpoints, max_ind, 4, 90, 110)

    bp_set = list(random.choice(bp_sets))

    lib = sr.Library.calc_from_config(bp_set, Decimal(1.0), lib_config)

    lib_json = lib.to_json()

    in_lib = sr.Library.from_json(lib_json)

    assert in_lib.breakpoints == lib.breakpoints
    assert in_lib.energy == lib.energy

    assert in_lib.energy_function.parents.alignment == \
        lib.energy_function.parents.alignment

    assert in_lib.mutation_rate == lib.mutation_rate
    assert in_lib.gg_prob == lib.gg_prob
    assert in_lib.gg_overhangs == lib.gg_overhangs

    # Should be fine since it was tested separately.
    assert in_lib.gg_enzyme.to_json() == lib.gg_enzyme.to_json()

    assert in_lib.amino_to_cdn == lib.amino_to_cdn


def test_dna_blocks(lib_config, breakpoints):
    parents = lib_config.energy_function.parents
    max_ind = len(parents.alignment)
    bp_sets = generate_bp_sets(breakpoints, max_ind, 4, 90, 110)

    bp_set = list(random.choice(bp_sets))

    lib = sr.Library.calc_from_config(bp_set, Decimal(1.0), lib_config)

    # TODO: Simulate Golden Gate and compare to parents.

    assert lib.dna_blocks
    for a in lib.dna_blocks:
        assert a
        for b in a:
            assert b


# TODO: Test find_best_overhangs.
