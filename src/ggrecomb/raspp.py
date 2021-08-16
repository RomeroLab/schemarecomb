"""Module for running the RASPP algorithm.

TODO: Refactor to separate codon_options, EnergyFunctions.
TODO: Add ABC for problem (RASPP superclass). Zheng et al. 2009 has joint
optimization of stability and diversity. Thus there are more algorithms that
may be useful to implement/replace RASPP. Does this necessitate renaming the
package to something like ggrecomb?
"""

from collections.abc import Iterable
from dataclasses import dataclass, field
from itertools import product, combinations
from typing import NamedTuple, Optional, Type, Union

from ggrecomb.energy_functions import EnergyFunction, SCHEMA
from ggrecomb.library import RecombinantLibrary
from ggrecomb import ParentAlignment
from ggrecomb import PDBStructure


'''
# List of optimal codons for E. coli.
AA_C31 = {'A': ['GCT', 'GCA'], 'R': ['CGT', 'CGA'], 'N': ['AAT'], 'D': ['GAT'],
          'C': ['TGT'], 'Q': ['CAA', 'CAG'], 'E': ['GAA'], 'G': ['GGT'],
          'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC'], 'L': ['TTA', 'TTG', 'CTA'],
          'K': ['AAA'], 'M': ['ATG'], 'F': ['TTT'], 'P': ['CCT', 'CCA'],
          'S': ['AGT', 'TCA'], 'T': ['ACA', 'ACT'], 'W': ['TGG'], 'Y': ['TAT'],
          'V': ['GTT', 'GTA'], '*': ['TGA'], '-': ['---']}
'''

# List of optimal codons for E. coli.
AA_C31 = {'A': ('GCT', 'GCA'), 'R': ('CGT', 'CGA'), 'N': ('AAT',),
          'D': ('GAT',), 'C': ('TGT',), 'Q': ('CAA', 'CAG'), 'E': ('GAA',),
          'G': ('GGT',), 'H': ('CAT', 'CAC'), 'I': ('ATT', 'ATC'),
          'L': ('TTA', 'TTG', 'CTA'), 'K': ('AAA',), 'M': ('ATG',),
          'F': ('TTT',), 'P': ('CCT', 'CCA'), 'S': ('AGT', 'TCA'),
          'T': ('ACA', 'ACT'), 'W': ('TGG',), 'Y': ('TAT',),
          'V': ('GTT', 'GTA'), '*': ('TGA',), '-': ('---',)}


class RASPP:
    """Recombination as a shortest path problem.

    Construct the graph depicted in Endelman et al 2004 figure 2. A node at
    index i in column c represents the selection of the (c-1)th breakpoint at
    the ith position in the multiple sequence alignment. Specifically, the ith
    amino acid in an aligned chimera will be the first residue of block c,
    while the (i-1)th amino acid is the last residue of block c-1. An edge
    between two nodes represents the block between the two breakpoints. A
    recombinant sequence library is therefore a traversal through the graph
    from node (c=0, i=0) to (n, i_f), where the sum of the edge weights in the
    path is the average energy of the library.

    Public attributes:
        columns: List of Lists of Nodes representing graph.
        breakpoints: Available recombination cross-over locations.

    Public methods:
        eval_bps: Average SCHEMA energy of sequence library with given
            breakpoints.
        eval_all_bps: eval_bps over all combinations of bps.
        min_bps: Find the set of breakpoints with minimum energy.
    """
    def __init__(
            self,
            parent_aln: ParentAlignment,
            n: int,
            start_overhang: Optional[tuple[int, str]] = None,
            end_overhang: Optional[tuple[int, str]] = None,
            codon_options: dict[str, tuple[str, ...]] = AA_C31,
            energy_func: Union[EnergyFunction, Type[EnergyFunction]] = SCHEMA
    ):
        """Init a RASPP graph.

        Args:
            parent_aln: Parent multiple sequence alignment.
            n: Number of breakpoints (crossovers) in output libraries.
            start_overhang: Positional shift and nucleotide sequence of Golden
                Gate site for vector insertion at start of sequence. Not
                factored into calculations if None.
            end_overhang: Positional shift and nucleotide sequence of Golden
                Gate site for vector insertion at end of sequence. Not factored
                into calculations if None.
            codon_options: Amino acid to available codon mapping. Used in
                Golden Gate site design and library sequence design. Change
                this to include or exclude certain codons based on codon
                optimization schemes, reassigned codons, etc.
            energy_func: Energy function used in calculations. See
                EnergyFunction for more information. May pass in class or
                instance.
        """

        # Calculate all crossover sites valid for Golden Gate Assembly.
        self.breakpoints, self._bp_to_group = _calculate_breakpoints(
            parent_aln, codon_options, start_overhang, end_overhang)

        if isinstance(energy_func, type):
            energy_func = energy_func(parent_aln)

        # Cache used for redundant <M> calculations. Certain adjacent
        # breakpoints will have the same <M>, e.g. if there is only one amino
        # acid in the alignment at each of the two breakpoint positions. Use
        # the self._bp_to_group map to convert a library's breakpoints to the
        # nonredundant breakpoint group indices, which will be the keys to
        # self._bp_group_M_cache.
        self._bp_group_M_cache: dict[tuple[int, ...], float] = {}

        self.pa = parent_aln

        # Build RASPP graph nodes.
        if start_overhang is None:
            start_overhang_list = None
        else:
            start_overhang_list = [start_overhang]
        first_col = [_Node(0, 0, start_overhang_list)]
        self.columns = [first_col]
        for col in range(1, n + 1):
            # No edges into node => exclude it.
            min_index = self.columns[col-1][0].index + 1
            new_col = [_Node(col, index, self.breakpoints[index])
                       for index in self.breakpoints if index >= min_index]
            self.columns.append(new_col)

        # Fill out edges between nodes.
        self.columns[0][0].fill_out(self.columns[1], energy_func)
        for col, next_col in zip(self.columns[1:], self.columns[2:]):
            for node in col[:-1]:
                node.fill_out(next_col, energy_func)

    # TODO: Refactor or remove next two methods to be compatible?
    def eval_bps(self, bps: Iterable[int]) -> float:
        """Average SCHEMA energy of sequence library with given breakpoints."""
        e = 0.0
        curr_node = self.columns[0][0]
        for bp in bps:
            edge = curr_node.out_edges[bp]
            e += edge.e_diff
            curr_node = edge.out_node
        return e

    def eval_all_bps(self) -> dict[tuple[int], float]:
        """eval_bps over all combinations of breakpoints."""
        valid_bps = list(self.breakpoints.keys())
        try:
            del valid_bps[0]
            del valid_bps[-1]
        except KeyError:
            pass
        bps_e = {bps: self.eval_bps(bps) for bps
                 in combinations(valid_bps, len(self.columns)-1)}
        return bps_e

    def min_bps(
        self,
        min_block_len: Optional[int] = None,
        max_block_len: Optional[int] = None
    ) -> 'RecombinantLibrary':
        """Find the set of breakpoints with minimum energy."""

        n = len(self.columns) - 1  # number of breakpoints
        N = len(str(self.pa.aligned_sequences[0].seq))

        if min_block_len is None:
            min_block_len = 0
        if max_block_len is None:
            max_block_len = N

        # Memoized traversal over columns of graph to find minimum library.
        # Record the minimum library for each node during traversal.
        first_node = self.columns[0][0]
        first_lib = _TempLibrary(0.0, first_node.breakpoint)
        next_col_min_libs = {first_node: first_lib}
        while next_col_min_libs:
            curr_col_min_libs = next_col_min_libs

            # print(list(curr_col_min_libs)[0].col)
            if list(curr_col_min_libs)[0].col == n:
                # We've reached the last column.
                break

            # Mapping from each node in the next column to it's minimum energy
            # library.
            next_col_min_libs = {}

            # Iterate over every node and it's minimum energy library in the
            # current column to find the minimum energy library for each node
            # in the next column.
            for curr_node, curr_min_lib in curr_col_min_libs.items():
                # Iterate over all edges from curr_node. Each edge results in
                # a candidate library built from curr_lib and the breakpoint
                # from the node at the end of the edge.
                for edge in curr_node.out_edges:
                    new_node = edge.out_node

                    # Check that the new block is of valid size.
                    next_block_len = new_node.index - curr_node.index
                    if not min_block_len <= next_block_len <= max_block_len:
                        continue

                    # Get the best library for new_node so far and compare to
                    # library from curr_node. If no best library or curr_node
                    # library is best, set as the best library.
                    best_lib = next_col_min_libs.get(new_node)
                    if best_lib is None or curr_min_lib.energy + edge.e_diff \
                       < best_lib.energy:
                        new_lib = curr_min_lib.expand(edge.e_diff,
                                                      new_node.breakpoint)
                        next_col_min_libs[new_node] = new_lib

        candidate_libs = [lib for node, lib in curr_col_min_libs.items()
                          if min_block_len <= N - node.index <= max_block_len]
        if not candidate_libs:
            raise LibrariesNotFoundError
        return_lib = min(candidate_libs, key=lambda lib: lib.energy)
        if not return_lib.breakpoints:
            raise LibrariesNotFoundError

        # add last breakpoint
        last_bp = {N: self.breakpoints[N]}
        temp_lib = return_lib.expand(0.0, last_bp)
        e = temp_lib.energy
        bps = temp_lib.breakpoints
        return_lib = RecombinantLibrary(e, bps, self.pa, self._bp_to_group,
                                        self._bp_group_M_cache)
        return return_lib

    def vary_m_proxy(self, L_min, L_max):
        """Choose E-optimal libraries by varying block length as an M proxy.

        TODO: Erase/refactor to get rid of old/naive min_bps() and
        calc_average_m().
        """
        chosen_libs = []

        # TODO: Can probably do problem reduction for each min_block_len
        # iteration. E.g. For a given min_block_len iteration, iterating
        # through max_block_len backwards allows us to say the most recently
        # found library is the minimum E until its max_block_len is greater
        # max_block_len.
        # TODO: DOUBLE CHECK THIS
        for min_block_len in range(L_min, L_max):
            curr_best = None
            for max_block_len in range(L_max, min_block_len, -1):
                print(min_block_len, max_block_len, '\r', end='')
                if curr_best is not None:
                    if curr_best.max_block_len <= max_block_len:
                        continue
                    else:
                        curr_best = None
                try:
                    lib = self.min_bps(min_block_len, max_block_len)
                    curr_best = lib

                    if lib not in chosen_libs:
                        lib.calc_average_m(self.pa, self._bp_group_map,
                                           self._bp_group_M_cache)
                        chosen_libs.append(lib)
                except LibrariesNotFoundError:
                    pass

        return chosen_libs

    def all_libraries(
        self,
        min_block_len: Optional[int] = None,
        max_block_len: Optional[int] = None
    ) -> list['RecombinantLibrary']:
        """Get all possible libraries with depth-first graph traversal."""
        libs = []
        n = len(self.columns) - 1  # number of breakpoints
        N = sorted(self.breakpoints)[-1]  # length of sequences

        if min_block_len is None:
            min_block_len = 0
        if max_block_len is None:
            max_block_len = N

        first_node = self.columns[0][0]
        first_lib = RecombinantLibrary(0.0, first_node.breakpoint)
        stack = [(first_node, first_lib)]
        while stack:
            curr_node, curr_lib = stack.pop()
            if curr_node.col == n:
                next_block_len = N - curr_node.index
                if min_block_len <= next_block_len <= max_block_len:
                    last_bp = {N: self.breakpoints[N]}
                    curr_lib = curr_lib.expand(0.0, last_bp)
                    curr_lib.calc_average_m(self.pa, self.bp_group_map)
                    libs.append(curr_lib)
                continue
            for edge in reversed(curr_node.out_edges.values()):
                next_block_len = edge.out_node.index - curr_node.index
                if min_block_len <= next_block_len <= max_block_len:
                    new_node = edge.out_node
                    new_lib = curr_lib.expand(edge.e_diff, new_node.breakpoint)
                    stack.append((new_node, new_lib))
        if not libs:
            msg = f'No libraries found with min_block_len={min_block_len}' \
                  f' and max_block_len={max_block_len}.'
            raise LibrariesNotFoundError(msg)
        return libs


@dataclass
class _TempLibrary:
    energy: float
    breakpoints: dict[int, list[tuple[int, str]]] = field(default_factory=dict)

    def expand(self, e_diff: float, new_bp: dict[int, list[tuple[int, str]]]):
        new_e = self.energy + e_diff
        new_bps = self.breakpoints | new_bp
        return type(self)(new_e, new_bps)


class LibrariesNotFoundError(Exception):
    pass


def _get_valid_patterns(
    cdn_sets: tuple[set[str]],
    reverse: bool = False
) -> tuple[set[str]]:
    """DNA patterns in the codon alignment for Golden Gate site construction.

    Helper function for _calculate_breakpoints. Each input set corresponds to
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

    TODO: Eliminate crossovers based on invalid/bad GG sites.
    """

    # If your codons are not three bases long, you're either wrong or doing
    # xenobiology.
    assert all(all(len(cdn) == 3 for cdn in cdn_set) for cdn_set in cdn_sets)

    # If not reverse, cdn_shrink removes the right side letter of each codon.
    # If reverse, cdn_shrink removes the left side letter of each codon.
    fwd_lim, bwd_lim = int(reverse), 1 - int(reverse)

    def cdn_shrink(cdns):
        return {cdn[fwd_lim:len(cdn)-bwd_lim] for cdn in cdns}

    # 1 AA per codon, so len 3 patterns exist iff len(cdn_sets) == 1.
    if len(cdn_sets) == 1:
        len_3_sets = tuple(cdn_sets)[0]
    else:
        len_3_sets = {}

    # Remove a codon letter and take intersection to get all len 2 patterns.
    cdn_sets = [cdn_shrink(cdns) for cdns in cdn_sets]
    len_2_sets = set.intersection(*cdn_sets)

    # Remove another letter and take intersection to get all len 1 patterns.
    cdn_sets = [cdn_shrink(cdns) for cdns in cdn_sets]
    len_1_sets = set.intersection(*cdn_sets)

    return len_1_sets, len_2_sets, len_3_sets


def _calculate_breakpoints(
    parent_aln: ParentAlignment,
    codon_options: dict[str, tuple[str, ...]],
    start_overhang: Optional[tuple[int, str]] = None,
    end_overhang: Optional[tuple[int, str]] = None,
) -> dict[int, list[tuple[int, str]]]:
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
    codon of the breakpoint (codon index b_k - 1).

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

    Args:
        parent_aln: parent alignment for breakpoint calculation.
        codon_options: Amino acid to available codon mapping. Used in
            Golden Gate site design and library sequence design. Change
            this to include or exclude certain codons based on codon
            optimization schemes, reassigned codons, etc.
        start_overhang: Positional shift and nucleotide sequence of Golden
            Gate site for vector insertion at start of sequence. Not
            factored into calculations if None.
        end_overhang: Positional shift and nucleotide sequence of Golden
            Gate site for vector insertion at end of sequence. Not factored
            into calculations if None.

    Returns:
        breakpoints: Mapping from breakpoint position to valid breakpoint
            overhangs.
        bp_to_group: Mapping from breakpoint index to index of nonredundant
            breakpoint groups.
    """

    # Check that codon_options contain all AAs in parent alignment.
    print('hello')
    print(parent_aln)
    aa_alphabet: list[tuple[str, ...]] = set().union(*parent_aln)
    print(aa_alphabet)
    assert all(aa in codon_options for aa in aa_alphabet)

    # Sets of codons for each amino acid at an alignment position.
    aln_cdns = [{codon_options[aa] for aa in pos_aas}
                for pos_aas in parent_aln]

    breakpoints = {}

    if start_overhang is not None:
        breakpoints[0] = [start_overhang]
        internal_bp_start = 1
    else:
        internal_bp_start = 0

    # Search possible Golden Gate sites by iterating over adjacent codons.
    for bp, (cdns1, cdns2) in enumerate(zip(aln_cdns, aln_cdns[1:]), 1):

        # Can't do Golden Gate at a site with gaps.
        # TODO: Do further Golden Gate validation.
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
            for pat1, pat2 in product(pat_set1, pat_set2):
                oh = (pos, pat1 + pat2)
                overhangs.append(oh)

        if not overhangs:
            # No valid overhangs, skip.
            continue

        breakpoints[bp] = overhangs

    if end_overhang is not None:
        breakpoints[len(aln_cdns)] = [end_overhang]
        internal_bp_end = len(breakpoints) - 1
    else:
        internal_bp_end = len(breakpoints)

    # TODO: Check if there's enough breakpoints for the run params.
    # Group the redundant breakpoints.
    breakpoint_groups = []

    # Iterate over the internal (not added start or end) breakpoints.
    bps = iter(sorted(breakpoints)[internal_bp_start:internal_bp_end])
    curr_bps = [next(bps)]
    for bp in bps:
        # If bp is not adjacent to last or has more than one amino acid, it
        # belongs to a new group.
        if bp != curr_bps[-1] + 1 or len(set(parent_aln[bp-1])) != 1:
            breakpoint_groups.append(curr_bps)
            curr_bps = [bp]
        else:
            curr_bps.append(bp)
    breakpoint_groups.append(curr_bps)

    # Get the mapping from each breakpoint index to the breakpoint's group
    # index.
    bp_to_group = {}
    for group_index, bp_group in enumerate(breakpoint_groups):
        for bp in bp_group:
            bp_to_group[bp] = group_index

    return breakpoints, bp_to_group


class _Edge(NamedTuple):
    """Edge between nodes in RASPP graph.

    Public Attributes:
        e_diff: Difference in energy between out_node and in_node.
        in_node: Input node.
        out_node: Output_node.
    """
    e_diff: float
    in_node: '_Node'
    out_node: '_Node'


@dataclass(repr=False)
class _Node:
    """Node in RASPP graph.

    Public Attributes:
        col: Column where this node is located.
        index: Position of breakpoint represented by this node.
        in_edges: Edges that point to this node.
        out_edges: Mapping from edge.index to edge for each edge with this node
            as in_node.

    Public Methods:
        fill_out: Fill out edges from this node to next column.
    """
    col: int
    index: int
    overhangs: Optional[list[tuple[int, str]]]
    in_edges: list[_Edge] = field(default_factory=list)
    out_edges: list[_Edge] = field(default_factory=list)

    @property
    def breakpoint(self):
        return {self.index: self.overhangs}

    def fill_out(
        self,
        next_col: list['_Node'],
        energy_func: EnergyFunction
    ) -> None:
        """Fill out edges from this node to next column.

        Args:
            next_col: Next column in graph.
            energy_func: Used to calculate energy between graph nodes.
        """

        def add_edge(energy, target_node):
            """Add an edge from self to target_node."""
            edge = _Edge(e, self, target_node)
            # self.out_edges[target_node.index] = edge
            self.out_edges.append(edge)
            target_node.in_edges.append(edge)

        # Next_col must have sorted indices.
        assert all(n1.index < n2.index for n1, n2
                   in zip(next_col, next_col[1:]))

        # Iterate on column, calculating energy and adding edges as needed.
        col_iter = iter(next_col)

        # Find a node such that the index is greater than this node's.
        # TODO: account for maxBL and minBL
        curr_node = next(col_iter)
        while curr_node.index <= self.index:
            curr_node = next(col_iter)

        # Get the energy of the first valid node in the column, make edge.
        e = energy_func.block_energy(self.index, curr_node.index)
        add_edge(e, curr_node)

        # Start iterating over breakpoint columns starting at next index.
        start_index = curr_node.index + 1

        # Get the next node. If this is a StopIteration, we are at the end of
        # the column, so we can return.
        try:
            curr_node = next(col_iter)
        except StopIteration:
            return

        for i in range(start_index, next_col[-1].index):
            # Calculate energy difference using previous energy and energy
            # difference between previous and current indicies.
            e += energy_func.increment_block_energy(self.index, i)

            # If we're at curr_node's index, we can add an edge to the node.
            if i == curr_node.index:
                add_edge(e, curr_node)
                curr_node = next(col_iter)

        # Add the last edge. Doing it here avoids the StopIteration.
        add_edge(e, curr_node)

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        in_edges = [e.in_node.index for e in self.in_edges]
        out_edges = [e.out_node.index for e in self.out_edges]
        ret = f'Node({self.col}, {self.index}, {self.breakpoint}, ' \
            f'{in_edges=}, {out_edges=})'
        return ret


if __name__ == '__main__':
    import sys

    loc = '../../tests/bgl3_sample/'
    pa = ParentAlignment.from_fasta(loc+'bgl3_sequences.fasta')
    pdb = PDBStructure.from_pdb_file(loc+'1GNX.pdb')
    pa.pdb_structure = pdb

    sys.exit()

    '''
    loc = '../tests/bgl3_sample/truncated/'
    pa = ParentAlignment.from_fasta(loc+'trunc.fasta')
    pdb = PDBStructure.from_pdb_file(loc+'trunc.pdb')
    '''
    loc = '../../tests/bgl3_sample/'
    pa = ParentAlignment.from_fasta(loc+'bgl3_sequences.fasta')
    pdb = PDBStructure.from_pdb_file(loc+'1GNX.pdb')
    pa.pdb_structure = pdb

    vector_overhangs = [(0, 'TATG'), (3, 'TGAG')]
    n = 4
    raspp = RASPP(pa, n, vector_overhangs[0], vector_overhangs[1])

    '''
    aln_sets = [set(a) for a in pa]
    breakpoint_groups = []
    bps = iter(sorted(raspp.breakpoints)[1:-1])
    curr_bps = [next(bps)]
    for bp in bps:
        print(bp, aln_sets[bp-1], aln_sets[bp])
        if bp != curr_bps[-1] + 1 or len(aln_sets[bp-1]) != 1:
            breakpoint_groups.append(curr_bps)
            curr_bps = [bp]
        else:
            curr_bps.append(bp)
    breakpoint_groups.append(curr_bps)
    print(breakpoint_groups)
    '''

    minBL, maxBL = 70, 105
    # minBL, maxBL = 30, 105
    import time

    # varied
    s = time.time()
    libs = raspp.vary_m_proxy(minBL, maxBL)
    '''
    print('\nvary_m_proxy')
    for rl in sorted(libs, key=lambda x: x.energy):
        print(rl.energy, list(rl.breakpoints), rl.min_block_len,
              rl.max_block_len, rl.average_m)
    '''
    print(time.time() - s)

    # all
    print('\nall_libraries')
    s = time.time()
    all_libs = raspp.all_libraries(minBL, maxBL)
    print('here')
    print(len(all_libs))
    for i, rl in enumerate(sorted(all_libs, key=lambda x: x.energy)):
        print(i, rl.energy, list(rl.breakpoints), rl.min_block_len,
              rl.max_block_len, rl.average_m)
    print(time.time() - s)
    1/0

    import matplotlib.pyplot as plt
    lib_me = set((rl.average_m, rl.energy) for rl in libs)
    alib_me = set((rl.average_m, rl.energy) for rl in all_libs)
    lib_plot = lib_me - alib_me
    alib_plot = alib_me - lib_me
    both_plot = lib_me & alib_me
    print(lib_plot)
    print(alib_plot)
    print(both_plot)

    all_breakpoints = [sorted(rl.breakpoints) for rl in all_libs]
    for rl in libs:
        if sorted(rl.breakpoints) not in all_breakpoints:
            print(sorted(rl.breakpoints))

    if lib_plot:
        plt.scatter(*zip(*lib_plot), label='varied')
    if alib_plot:
        plt.scatter(*zip(*alib_plot), label='all')
    plt.scatter(*zip(*both_plot), label='both')
    plt.xlabel('<M>')
    plt.ylabel('<E>')
    plt.legend()
    plt.show()

    """
    rl = min(all_libs, key=lambda x: x.energy)
    print(rl)
    print(rl.min_block_len, rl.max_block_len)
    """

    '''
    for minBL, maxBL in [(1, 8), (2, 5), (None, None)]:
        all_libs = raspp.all_libraries(minBL, maxBL)
        print(minBL, maxBL, [sorted(lib.breakpoints) for lib in all_libs])
    '''
    '''
    raspp.eval_all_bps()
    raspp.min_bps()

    N = len(pa.aligned_sequences[0].seq)
    opt_time = 0.0
    opt_single_time = 0.0
    naive_time = 0.0
    for bps in combinations(range(1, 9), 5):
        m = average_m(pa, (0,) + bps)
        print(bps, m)
    '''
