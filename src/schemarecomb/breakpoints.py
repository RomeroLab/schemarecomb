"""Recombination site discovery in :class:`~schemarecomb.ParentSequences`."""

import itertools
from typing import NamedTuple
from typing import Optional

from schemarecomb import ParentSequences


class Overhang(NamedTuple):
    """Valid Golden Gate Assembly overhangs (sticky ends).

    Parameters:
        ind: Index of the start site relative to the codon at one less than the
            breakpoint's position - 1. For example, the (1, 'TGTG') overhang
            would be a sticky end that is one base pair right of the codon at
            position-1. So if the alignment only has 'M' and 'W' at position-1
            and position, respectively, "ATG" and "TGG" are the only codons
            that may appear at these indices. Then, with the overhang as
            upper-case, the assembled DNA sequence will be "...aTGTGg...". An
            ind of -1 implies that the overhang is a dummy.
        seq: Sequence of the sticky end.

    """
    ind: int
    seq: str


class BreakPoint(NamedTuple):
    """Valid site for Golden Gate Assembly in a parent alignment.

    A BreakPoint is where an alignment may be cut and recombined to form
    chimeric sequences. Also known as "crossover" or "recombination site".

    Parameters:
       position: Index within a ParentSequence.alignment where the BreakPoint
            is located. Note carefully that this attribute denotes where a new
            chimeric block begins. The last amino acid of the previous block
            is located at position-1. The Golden Gate sticky ends will overlap
            the codons of amino acids position-1 and position.
        overhangs: Valid Golden Gate overhangs, containing the overhang DNA
            start site relative to the codon at position-1 and the overhang
            string. For example, the (1, 'TGTG') overhang would be a sticky end
            that is one base pair right of the codon at position-1. So if the
            alignment only has 'M' and 'W' at position-1 and position,
            respectively, "ATG" and "TGG" are the only codons that may appear
            at these indices. Then, with the overhang as upper-case, the
            assembled DNA sequence will be "...aTGTGg...".

    """
    position: int
    overhangs: list[Overhang]


def _get_valid_patterns(
    cdn_sets: list[set[str]],
    reverse: bool = False
) -> tuple[set[str], set[str], set[str]]:
    """DNA patterns in the codon alignment for Golden Gate site construction.

    Helper function for calculate_breakpoints. Each input set corresponds to
    the valid codons that encode a given amino acid at a given position in the
    parent alignment.  For example, if a certain position has isoleucine and
    threonine, and the valid codons for these amino acids are ('ATT', 'ATC')
    and ('ACC', 'ACT'), respectively, ({'ATT', 'ATC'), {'ACC', 'ACT'}) will be
    passed in. The length 1, 2, and 3 patterns are returned. See the example
    below for more information.

    Args:
        cnd_sets: Collections of codons that code for each amino acids at a
            given site in the parent alignment.
        reverse: If the patterns should be explored right to left instead of
            right to left (default).

    Returns:
        Collections of all length 1, 2, and 3 patterns found.

    Example:
        _get_valid_patterns(({'ATT', 'ATC'}, {'ACC', 'ACT'})) will return
            ({'A'}, {}, {}) since 'A' can be found as the first letter at least
            once in each set. There are no length 2 patterns because there are
            no common letters at the second codon position.
        _get_valid_patterns(({'ATT', 'ATC'}, {'ACC', 'ACT'}), True) will return
            ({'C', 'T'}, {}, {}) since 'C' and 'T' can both be found as the
            last letter at least once in each set.
    """

    if not cdn_sets:
        raise ValueError('cdn_sets must not be empty.')

    # If your codons are not three bases long, you're either wrong or doing
    # xenobiology.
    if not all(all(len(cdn) == 3 for cdn in cdn_set) for cdn_set in cdn_sets):
        raise ValueError('Codons should be three bases long.')

    # If not reverse, cdn_shrink removes the right side letter of each codon.
    # If reverse, cdn_shrink removes the left side letter of each codon.
    fwd_lim, bwd_lim = int(reverse), 1 - int(reverse)

    def cdn_shrink(cdns: set[str]) -> set[str]:
        return {cdn[fwd_lim:len(cdn)-bwd_lim] for cdn in cdns}

    # 1 AA per codon, so len 3 patterns exist iff len(cdn_sets) == 1.
    if len(cdn_sets) == 1:
        return set(), set(), set(cdn_sets[0])

    # Remove a codon letter and take intersection to get all len 2 patterns.
    cdn_sets = [cdn_shrink(cdns) for cdns in cdn_sets]
    len_2_sets = set.intersection(*cdn_sets)

    # Remove another letter and take intersection to get all len 1 patterns.
    cdn_sets = [cdn_shrink(cdns) for cdns in cdn_sets]
    len_1_sets = set.intersection(*cdn_sets)

    return len_1_sets, len_2_sets, set()


def calculate_breakpoints(
    parents: ParentSequences,
    codon_options: dict[str, set[str]],
    start_overhangs: Optional[list[Overhang]] = None,
    end_overhangs: Optional[list[Overhang]] = None,
) -> list[BreakPoint]:
    """Calculate the breakpoints for the parent alignment.

    Currently supports length four Golden Gate sites.

    # TODO: move this to module header?
    A breakpoint is potential site for recombination of the parental sequences.
    The terms "breakpoint" and "crossover" are used synonymously. A
    breakpoint's position is the index of the breakpoint's second amino acid in
    the parent alignment. This definition allows the breakpoint position to
    nicely slice the parent alignment into blocks. For example, let the
    position of the k-1th, kth, and k+1th breakpoints be b_k_1, b_k, b_k__1,
    respectively. Then for parent sequence p, the block before the breakpoint
    is p[b_k_1:b_k] and the block after is p[b_k:b_k__1]. Valid breakpoints
    contain candidate overhangs, defined below.

    Overhangs are the DNA sticky ends used in a Golden Gate reaction. In this
    module, each overhang consists of a positional shift and DNA sequence. The
    former is overhang's position relative to the first base pair in the left
    codon of the breakpoint (codon index b_k - 1). See :class:`BreakPoint`.

    Example: Suppose the amino acids at indices b_k-1 and b_k at 'I' and 'T'
        for all parents in the parent alignment and codon_options contains
        {'M': ('ATG,), 'T': ('ACC')}}. Then for each parent p, the only valid
        codon sequence for p[b_k-1:b_k+1] is ATGACC, and the overhangs are
        (0, 'ATGA'), (1, 'TGAC'), (2, 'GACC'). This breakpoint will be
        represented in the breakpoints dictionary as
        {b_k: [(0, 'ATGA'), (1, 'TGAC'), (2, 'GACC')]}.

    Certain adjacent breakpoints will result in the same <M>, e.g. if alignment
    positions i and i+1 only consist of {'A'} and {'C'} respectively, the
    libraries of breakpoints (..., i, ...) and (..., i+1, ...) will have the
    same <M> if all other breakpoints are the same. Therefore, this function
    groups the redundant breakpoints and returns a mapping from the original
    breakpoint index to the nonredundant group index.

    Parameters:
        parents: parent alignment for breakpoint calculation.
        codon_options: Amino acid to available codon mapping. Used in
            Golden Gate site design and library sequence design. Change
            this to include or exclude certain codons based on codon
            optimization schemes, reassigned codons, etc.
        start_overhangs: Positional shift and nucleotide sequence of Golden
            Gate site for vector insertion at start of sequence. Not
            factored into calculations if None.
        end_overhangs: Positional shift and nucleotide sequence of Golden
            Gate site for vector insertion at end of sequence. Not factored
            into calculations if None.

    Returns:
        Mapping from breakpoint position to valid breakpoint overhangs.

    """

    try:
        alignment = parents.alignment
    except AttributeError:
        raise ValueError('Input ParentSequences must be aligned to calculate '
                         'valid breakpoints.')

    # Check that codon_options contain all AAs in parent alignment.
    aa_alphabet = set(itertools.chain(*alignment))
    assert all(aa in codon_options for aa in aa_alphabet)

    # Sets of codons for each amino acid at an alignment position.
    aln_cdns = [[codon_options[aa] for aa in pos_aas] for pos_aas in alignment]

    breakpoints = []

    if start_overhangs is not None:
        start_bp = BreakPoint(0, start_overhangs)
        breakpoints.append(start_bp)

    # Search possible Golden Gate sites by iterating over adjacent codons.
    for bp_index, (cdns1, cdns2) in enumerate(zip(aln_cdns, aln_cdns[1:]), 1):

        # Can't do Golden Gate at a site with gaps.
        if '---' in cdns1 or '---' in cdns2:
            continue

        # Get the [len 1], [len 2], and [len 3] patterns for the two sites
        # as list of lists.
        pos1_bins = _get_valid_patterns(cdns1, reverse=True)
        pos2_bins = _get_valid_patterns(cdns2)

        overhangs = []

        # Combine the pattern sets to get all the len 4 overhangs. The
        # pos1 list is reversed so that the len 3 patterns go with the
        # pos2 len 1 patterns, len 2 with len 2, etc.
        pattern_iter = zip(reversed(pos1_bins), pos2_bins)
        for pos, (pat_set1, pat_set2) in enumerate(pattern_iter):
            for pat1, pat2 in itertools.product(pat_set1, pat_set2):
                oh = Overhang(pos, pat1 + pat2)
                overhangs.append(oh)

        if not overhangs:
            # No valid overhangs, skip.
            continue

        bp = BreakPoint(bp_index, overhangs)
        breakpoints.append(bp)

    if end_overhangs is not None:
        end_bp = BreakPoint(len(aln_cdns), end_overhangs)
        breakpoints.append(end_bp)

    return breakpoints


def block_indices(
    breakpoints: list[BreakPoint],
    aln_len: int
) -> list[tuple[int, int]]:
    """Get the blocks' start and end indices given a list of breakpoints.

    The start index is inclusive and the end index is exclusive, i.e. slice
    indexing. For example, if a block starts at 40 and ends at 83, a parent
    sequence can be sliced by these indices to result in a chimeric block as
    such: p2_block1 = parent2[40:83]. This also implies that the end index of
    this block is the start index of the next block, e.g. p2_block2 =
    parent2[83:132].

    This function also handles the presence or absence of breakpoints at the
    zero position and alignment length position. These positions usually have
    fixed breakpoints (across all libraries) corresponding to the target vector
    that will hold the assembled chimera, but may be omitted depending on
    application. Regardless, this function will temporarily add breakpoint
    indices at the first and last positions (0 and aln_len) if necessary to
    correctly compute the chimeric blocks' start and end indicies.

    Parameters:
        breakpoints: BreakPoints where the parent alignment is cut, usually
            comes from a library.
        aln_len: The length of the alignment. This is the end index of the
            last block, which may already be in breakpoints. Generally passed
            in as len(parents.alignment), where parents is a ParentSequences.

    Returns:
        The (start, end) indices of each chimeric block.

    Raises:
        ValueError: If breakpoints is empty.

    """
    if not breakpoints:
        raise ValueError('breakpoints must not be empty.')

    # handle if breakpoints is missing 0 index or aln_len index.
    bp_inds = tuple(sorted(bp.position for bp in breakpoints))
    if bp_inds[0] != 0:
        bp_inds = (0,) + bp_inds
    if bp_inds[-1] != aln_len:
        bp_inds += (aln_len,)

    return [(i, f) for i, f in zip(bp_inds, bp_inds[1:])]
