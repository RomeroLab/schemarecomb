import pytest

import schemarecomb as sr
from schemarecomb import breakpoints


def test_get_valid_patterns():
    ret = breakpoints._get_valid_patterns(({'ATT', 'ATC'}, {'ACC', 'ACT'}))
    assert ret == ({'A'}, set(), set())

    ret = breakpoints._get_valid_patterns(({'ATT', 'ATC'}, {'ACC', 'ACT'}),
                                          True)
    assert ret == ({'C', 'T'}, set(), set())

    with pytest.raises(ValueError):
        breakpoints._get_valid_patterns(({'ACT', 'ACG'}, {'TCA', 'TA'}))

    with pytest.raises(ValueError):
        breakpoints._get_valid_patterns(())

    ret = breakpoints._get_valid_patterns(({'ATT', 'ATC'},))
    assert ret == (set(), set(), {'ATT', 'ATC'})


def test_calculate_breakpoints(bgl3_parents_aln, AA_C31):
    start = [(2, 'TATG')]
    end = [(0, 'TGAG')]
    sr.breakpoints.calculate_breakpoints(bgl3_parents_aln, AA_C31, start, end)
