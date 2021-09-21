"""Restriction enzymes for Golden Gate Assembly.

This module handles most of the Golden Gate Assembly computation. gg_prob is an
estimation of Golden Gate reaction efficiency, which may vary depending on the
selection of overhangs that must anneal. The biochemical theory on which the
module is based on was developed independently, but the data used and a
more detailed explanation can be found in `Pryor, Potapov, et al. 2020
<https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0238592>`_.

"""

from collections import deque
from decimal import Decimal
import itertools
import json
from math import prod
from operator import attrgetter
import os
from random import choices
from typing import NamedTuple

from schemarecomb.breakpoints import BreakPoint
from schemarecomb.breakpoints import Overhang


def read_ligation_data(fn: str) -> dict[tuple[str, str], Decimal]:
    """Read in ligation count data.

    The file must have annealing counts for each pair of overhangs. See the
    $PKG_ROOT/src/schemarecomb/gg_data directory for examples. The data in this
    directory comes from `Pryor, Potapov, et al. 2020 <https://journals.plos.
    org/plosone/article?id=10.1371/journal.pone.0238592>`_, but directly
    converted to CSV for ease of input.

    Parameters:
        fn: Name of data file.

    Returns:
        Mapping from two overhang strings to number of ligations observed
            between the two overhangs.

    """
    with open(fn) as f:
        headers = next(f).strip().split(',')
        lines = [line.strip().split(',') for line in f]

    overhangs = headers[1:]
    golden_gate_counts = {}

    for overhang1, *gg_data in lines:
        for overhang2, ligation_count in zip(overhangs, gg_data):
            count = Decimal(ligation_count)
            golden_gate_counts[overhang1, overhang2] = count

    return golden_gate_counts


_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def _reverse_complement(dna_str):
    return ''.join(reversed([_complement[bp] for bp in dna_str]))


class _StackElement(NamedTuple):
    """Element of the max_gg_prob stack.

    Attributes:
        prob: gg_prob of current element.

    """
    prob: Decimal
    ohs: list[Overhang]
    oh_denoms: list[Decimal]


class RestrictionEnzyme:
    """Type IIS restriction enzyme used in Golden Gate Assembly reaction.

    Only restriction enzymes that create three or four base overhangs are
    currently supported.

    Parameters:
        name: Name of the restriction enzyme.
        recognition_seq: Sequence of DNA bases recognized by enzyme. By
            convention, the cut site must be downstream of this sequence.
        top_cut_dist: Distance from the end of recognition_seq where the top
            strand (containing the recognition_seq) is cut.
        bottom_cut_dist: Distance from the end of recognition_seq where the
            bottom (complementary) strand is cut.
        ligation_probs: For each pair of possible DNA overhangs, the relative
            count of ligating together. Generally comes from the
            read_ligation_data function.

    Attributes:
        name: Name of the restriction enzyme.
        recognition_seq (str): Sequence of DNA bases recognized by enzyme. By
            convention, the cut site must be downstream of this sequence.
        top_cut_dist (int): Distance from the end of recognition_seq where the
            top strand (containing the recognition_seq) is cut.
        bottom_cut_dist (int): Distance from the end of recognition_seq where
            the bottom (complementary) strand is cut.
        overhang_len (int): Length of the overhang (sticky end) sequence.
        ligation_counts: For each pair of possible DNA overhangs, the relative
            count of ligating together.

    Raises:
        ValueError: If the input overhang length is not 3 or 4, or if the keys
            of ligation_counts does not contain all possible overhangs of the
            correct length.

    """
    def __init__(
        self,
        name: str,
        recognition_seq: str,
        top_cut_dist: int,
        bottom_cut_dist: int,
        ligation_counts: dict[tuple[str, str], Decimal],
    ):
        # Only overhang length of 3 or 4 are supported.
        oh_len = abs(top_cut_dist - bottom_cut_dist)
        if abs(oh_len) not in (3, 4):
            raise ValueError('Overhang length must be 3 or 4.')

        # Check that all overhangs are present as keys in ligation_counts.
        all_ohs = {''.join(oh) for oh in itertools.product('ACGT',
                                                           repeat=oh_len)}
        valid_keys = set(itertools.product(all_ohs, repeat=2))
        if valid_keys != ligation_counts.keys():
            raise ValueError('ligation_counts must have all pairs of possible '
                             'overhangs as keys.')

        self.name = name
        self.recognition_seq = recognition_seq
        self.top_cut_dist = top_cut_dist
        self.bottom_cut_dist = bottom_cut_dist
        self.ligation_counts = ligation_counts

    @property
    def overhang_len(self):
        return abs(self.top_cut_dist - self.bottom_cut_dist)

    # TODO: document this
    @property
    def for_restriction_site(self):
        num_sep_bases = min(self.top_cut_dist, self.bottom_cut_dist)
        sep_bases = ''.join(choices('acgt', k=num_sep_bases))
        return self.recognition_seq.lower() + sep_bases

    @property
    def rev_restriction_site(self):
        num_sep_bases = min(self.top_cut_dist, self.bottom_cut_dist)
        sep_bases = ''.join(choices('acgt', k=num_sep_bases))
        rev_rec_site = _reverse_complement(self.recognition_seq)
        return sep_bases + rev_rec_site.lower()

    @classmethod
    def from_name(cls, name: str) -> 'RestrictionEnzyme':
        """Construct specific restriction enzyme from allowed names.

        The allowed values for name are 'BsaI-HFv2', 'BbsI-HF', 'BsmBI-v2',
        'Esp3I', 'SapI'. Search on the NEB website for more information about
        these enzymes.

        """
        enzyme_params = {
            'BsaI-HFv2': ('GGTCTC', 1, 5),
            'BbsI-HF': ('GAAGAC', 2, 6),
            'BsmBI-v2': ('CGTCTC', 1, 5),
            'Esp3I': ('CGTCTC', 1, 5),
            'SapI': ('GCTCTTC', 1, 4),
        }
        dir_path = os.path.dirname(os.path.realpath(__file__))
        fn = os.path.join(dir_path, 'gg_data', name + '.csv')
        params = enzyme_params[name]
        lig_counts = read_ligation_data(fn)
        return cls(name, *params, lig_counts)

    def lc(self, oh_seq1: str, oh_seq2: str) -> Decimal:
        """Ligation counts between oh1 and oh2."""
        return self.ligation_counts[(oh_seq1, oh_seq2)]

    def _init_stack_element(self, oh: Overhang) -> _StackElement:
        """Init method for the first set of overhangs in max_gg_prob."""
        oh_seq = oh.seq
        coh_seq = _reverse_complement(oh_seq)  # Complementary overhang.

        # Get three types of ligation counts.
        oh_oh = self.lc(oh_seq, oh_seq)
        oh_coh = self.lc(oh_seq, coh_seq)
        coh_coh = self.lc(coh_seq, coh_seq)

        # gg_prob of overhang and its complement.
        prob = oh_coh**2 / ((oh_coh+oh_oh)*(oh_coh+coh_coh))

        ohs = [oh, Overhang(-1, coh_seq)]
        oh_denoms = [oh_oh + oh_coh, coh_coh + oh_coh]

        return _StackElement(prob, ohs, oh_denoms)

    def gg_prob(self, overhangs: list[Overhang],
                has_complements: bool = False) -> Decimal:
        """Golden Gate Assembly efficiency."""
        oh_seqs = [oh.seq for oh in overhangs]
        if not has_complements:
            complements = [_reverse_complement(oh_seq) for oh_seq in oh_seqs]
            # overhangs += list(reversed(complements))
            # overhangs += complements
            new_oh_seqs = []
            for oh_seq, coh_seq in zip(oh_seqs, complements):
                new_oh_seqs.extend([oh_seq, coh_seq])
            oh_seqs = new_oh_seqs

        prob = Decimal(1.0)

        for oh1_seq in oh_seqs:
            tot = sum(self.ligation_counts[(oh1_seq, oh2_seq)] for oh2_seq
                      in oh_seqs)
            oh1_complement = _reverse_complement(oh1_seq)
            p_oh1 = self.ligation_counts[(oh1_seq, oh1_complement)] / tot
            prob *= p_oh1

        return prob

    def max_gg_prob(
        self,
        breakpoints: list[BreakPoint],
        thresh: Decimal = Decimal(1.0)
    ) -> tuple[Decimal, list[Overhang]]:
        r"""Calculate the optimal set of overhangs.

        For each BreakPoint, find the best overhang such that the overhang set
        has the largest gg_prob. The thresh parameter enables a shortcut for
        faster runtime. Example usage: when constructing many libraries, we
        might not care about the specific gg_prob, as long as each library can
        be shown to have a gg_prob that's good enough. In this case, you can
        set thresh to 0.95, as this indicates a 95% Golden Gate efficiency,
        which is plenty if you transform the assembly. This will speed up the
        most computationally intensive part of the process. Then when you've
        selected a library, you can call this method again with thresh=1.0.

        This method optimizes over the Cartesian product of the overhang sets
        in breakpoints. Although it's significantly faster than calling gg_prob
        while iterating over itertools.product, it's most likely still O(m^n),
        where m is the number of breakpoints in a chimera and n is the number
        of overhangs for each breakpoint. Therefore, for large problems (>5
        parents, >10 breakpoints) you probably want to use a heuristic, such as
        simulated annealing, which may or may not be implemented in future
        versions of this module.

        Implementation details:
            Say we want to find the gg_prob of the set of overhangs
            :math:`\{x_1,...,x_{n-1},x_n\}`. let :math:`p_{n-1}` be the gg_prob
            of the set of overhangs :math:`\{x_1,...,x_{n-1}\}` and
            :math:`L_{ij}` be the number of ligations between overhangs
            :math:`x_i` and :math:`x_j` in the Pryor, Potapov, et al. 2020
            data. Then

            .. math::

                p_n = \prod_{i=1}^n \dfrac{L_{ii}}{\sum_{j=1}^n L_{ij}}.

            This is equation 1 in Pryor, Potapov, et al. 2020. Let

            .. math::

                D_{ni} = \sum_{j=1}^{n-1} L_{ij}.

            Then

            .. math::

                p_n = p_{n-1} \prod_{i=1}^{n-1} \dfrac{D_{ni}}{D_{ni} + L_{in}}
                \dfrac{L_{nn}}{\sum_{j=1}^n L_{nj}}.

            So for each step :math:`n`, if we calculate and store :math:`p_n`
            and :math:`D_{ni}`, we can recursively optimize over the sets of
            overhangs. Since :math:`p_n \leq p_{n-1}`, we don't need to
            consider a given :math:`x_n` if :math:`p_n` is smaller than the
            best :math:`p_m` we've found so far, where :math:`m` is the number
            of total overhang sets.

            The need to consider the overhangs' complements complicates the
            actual implementation, but this is rectified through relatively
            straightforward changes to the equations above. Since the gg_prob
            for the next overhang and its complement can be calculated
            together, they can both be incorporated in a single recusive step.

        Parameters:
            breakpoints: Mapping from breakpoint index to BreakPoint. The
                Cartesian product of the overhang attribute of each BreakPoint
                will be optimized over.
            thresh: Minimum gg_prob for a complete set of overhangs to be
                considered optimal. This set will immediately be returned,
                even if a set with larger gg_prob may exist.

        Returns:
            - gg_prob of returned overhang set.
            - best overhang set found, or the first set found with gg_prob
              greater than thresh.

        Raises:
            ValueError: If len(breakpoints) < 2 or thresh is either more than
                1.0 or less than or equal to 0.0.
        """

        if len(breakpoints) < 2:
            raise ValueError('Length of breakpoints must be greater than 1.')

        if not 0.0 < thresh <= 1.0:
            raise ValueError('thresh must be more than 0.0 and no more than '
                             '1.0.')

        '''
        overhang_sets: list[list[tuple[int, str]]] = []
        for breakp in sorted(breakpoints, key=lambda bp: bp.position):
            oh_set = list(set(breakp.overhangs))
            overhang_sets.append(oh_set)
        '''

        sorted_bps = sorted(breakpoints, key=attrgetter('position'))
        overhang_sets = [bp.overhangs for bp in sorted_bps]

        max_gg_prob = Decimal(0.0)
        max_overhangs: list[Overhang]

        initial_set = overhang_sets[0]
        comp_stack = deque(self._init_stack_element(oh) for oh in initial_set)

        while comp_stack:
            prob, ohs, oh_denoms = comp_stack.pop()

            if prob <= max_gg_prob:
                continue

            stack_add = []

            # Get the index of the next overhang set.
            oh_index = len(ohs) // 2
            for next_oh in overhang_sets[oh_index]:
                noh_str = next_oh.seq
                ncoh_str = _reverse_complement(noh_str)

                # Next denominators to track. Not including denoms for next_oh
                # and next_coh.
                next_denoms = [d + self.lc(oh.seq, noh_str)
                               + self.lc(oh.seq, ncoh_str)
                               for oh, d in zip(ohs, oh_denoms)]

                denom_zip = zip(oh_denoms, next_denoms)
                D_prod = Decimal(prod(d / next_d for d, next_d in denom_zip))

                next_coh = Overhang(-1, ncoh_str)
                next_overhangs = ohs + [next_oh, next_coh]

                # Denominators for next_oh and next_coh.
                n_denom = sum(self.lc(oh.seq, noh_str) for oh
                              in next_overhangs)
                n_denom = Decimal(n_denom)  # avoids mypy error
                nc_denom = sum(self.lc(oh.seq, ncoh_str) for oh
                               in next_overhangs)
                nc_denom = Decimal(nc_denom)  # avoids mypy error
                next_denoms += [n_denom, nc_denom]

                noh_ncoh = self.lc(noh_str, ncoh_str)
                next_prob = prob * D_prod * noh_ncoh**2 / n_denom / nc_denom

                # next_prob < prob, so don't add to stack if we've already
                # found a better one.
                if next_prob > max_gg_prob:
                    if oh_index+1 == len(breakpoints):
                        # We're at the end of breakpoints, this is a valid set.
                        max_gg_prob = next_prob
                        max_overhangs = next_overhangs
                        if max_gg_prob >= thresh:
                            # Found sequences that are good enough, stop
                            # searching. Get rid of complements.
                            max_overhangs = max_overhangs[::2]
                            return max_gg_prob, max_overhangs
                    else:
                        new_elem = _StackElement(next_prob, next_overhangs,
                                                 next_denoms)
                        stack_add.append(new_elem)

            if stack_add:
                # Greedy, best element is top of stack.
                stack_add.sort(key=lambda e: e.prob)
                comp_stack.extend(stack_add)

        try:
            # Get rid of overhangs.
            max_overhangs = max_overhangs[::2]
            return max_gg_prob, max_overhangs
        except NameError:
            # max_overhangs does not exist.
            raise ValueError('Maximum overhang set could not be found.')

    def to_json(self) -> str:
        lig_counts = [list(k) + [str(v)] for k, v
                      in self.ligation_counts.items()]
        out_list = [
            self.name,
            self.recognition_seq,
            self.top_cut_dist,
            self.bottom_cut_dist,
            lig_counts,
        ]
        return json.dumps(out_list)

    @classmethod
    def from_json(cls, in_json: str):
        name, req_seq, top_dist, bot_dist, lig_counts = json.loads(in_json)
        lig_counts = {(k1, k2): Decimal(v) for k1, k2, v in lig_counts}
        return cls(name, req_seq, top_dist, bot_dist, lig_counts)
