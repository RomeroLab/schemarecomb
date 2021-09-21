from decimal import Decimal
from itertools import product
import random

import pytest

import schemarecomb as sr


@pytest.fixture
def breakpoints(bgl3_parents_aln, AA_C31):
    start = [sr.breakpoints.Overhang(2, 'TATG')]
    end = [sr.breakpoints.Overhang(0, 'TGAG')]
    return sr.breakpoints.calculate_breakpoints(bgl3_parents_aln, AA_C31,
                                                start, end)


def test_max_gg_prob(breakpoints):
    resenz = sr.restriction_enzymes.RestrictionEnzyme.from_name('BsaI-HFv2')

    start_bp = breakpoints[0]
    end_bp = breakpoints[-1]
    internal_bps = breakpoints[1:-1]

    for _ in range(1000):
        opt_bps = [start_bp] + random.sample(internal_bps, 5) + [end_bp]
        gg_prob, opt_ohs = resenz.max_gg_prob(opt_bps)

        # Ideal from max_gg_prob agrees with gg_prob for non-trivial problems.
        assert abs(gg_prob - resenz.gg_prob(opt_ohs)) < Decimal(10**-6)

        # The best found is among the best found by brute-force.
        oh_sets = [[oh for oh in bp.overhangs] for bp in opt_bps]
        curr_best = Decimal(0.0)
        for cand_ohs in product(*oh_sets):
            cand_prob = resenz.gg_prob(cand_ohs)
            if cand_prob > curr_best:
                curr_best = cand_prob
                brute_opts = {frozenset(cand_ohs)}
            elif cand_prob == curr_best:
                brute_opts.add(frozenset(cand_ohs))
        try:
            # Might find ideal in different order.
            assert set(opt_ohs) in brute_opts
        except NameError:
            raise ValueError('No valid candidate found.')


def test_from_name():
    for name in ('BsaI-HFv2', 'BbsI-HF', 'BsmBI-v2', 'Esp3I', 'SapI'):
        sr.restriction_enzymes.RestrictionEnzyme.from_name(name)


def test_json():
    resenz = sr.restriction_enzymes.RestrictionEnzyme.from_name('BsaI-HFv2')
    resenz_json = resenz.to_json()
    resenz_in = sr.restriction_enzymes.RestrictionEnzyme.from_json(resenz_json)
    assert resenz.name == resenz_in.name
    assert resenz.recognition_seq == resenz_in.recognition_seq
    assert resenz.top_cut_dist == resenz_in.top_cut_dist
    assert resenz.bottom_cut_dist == resenz_in.bottom_cut_dist
    assert resenz.ligation_counts == resenz_in.ligation_counts
