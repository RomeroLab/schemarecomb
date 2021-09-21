"""Module for running the RASPP algorithm."""

from dataclasses import dataclass
from dataclasses import field
from decimal import Decimal
from itertools import combinations
from operator import attrgetter
from typing import Iterable
from typing import NamedTuple
from typing import Optional
from typing import Type
from typing import Union

from schemarecomb import Library
from schemarecomb import ParentSequences
from schemarecomb.breakpoints import BreakPoint
from schemarecomb.breakpoints import calculate_breakpoints
from schemarecomb.breakpoints import Overhang
from schemarecomb.energy_functions import EnergyFunction
from schemarecomb.energy_functions import SCHEMA
from schemarecomb.libraries import LibraryConfig
from schemarecomb.libraries import MutationRateCache
from schemarecomb.restriction_enzymes import RestrictionEnzyme


@dataclass
class _TempLibrary:
    """Simpler representation of Library for RASPP purposes."""
    energy: Decimal
    breakpoints: list[BreakPoint] = field(default_factory=list)

    def expand(self, e_diff: Decimal, new_bp: BreakPoint):
        new_e = self.energy + e_diff
        new_bps = self.breakpoints + [new_bp]
        return type(self)(new_e, new_bps)

    def to_library(self, lib_config):
        # Need to remove first library if it doesn't have overhangs.
        bps = sorted(self.breakpoints, key=attrgetter('position'))
        if bps[0].position == 0 and bps[0].overhangs == []:
            bps = bps[1:]
        return Library.calc_from_config(bps, self.energy, lib_config)


class _Edge(NamedTuple):
    """Edge between nodes in RASPP graph.

    Public Attributes:
        e_diff: Difference in energy between out_node and in_node.
        in_node: Input node.
        out_node: Output_node.
    """
    e_diff: Decimal
    in_node: '_Node'
    out_node: '_Node'


@dataclass(repr=False)
class _Node:
    """Node in RASPP graph. Represents a BreakPoint.

    Public Attributes:
        col: Column where this node is located.
        position: Position of breakpoint represented by this node.
        bp: BreakPoint held by node.
        in_edges: Edges that point to this node.
        out_edges: Mapping from edge.position to edge for each edge with this
            node as in_node.
        breakpoint: BreakPoint represented by this node.

    Public Methods:
        fill_out: Fill out edges from this node to next column.

    """
    col: int
    bp: BreakPoint
    in_edges: list[_Edge] = field(default_factory=list)
    out_edges: list[_Edge] = field(default_factory=list)

    @property
    def position(self):
        return self.bp.position

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
            # self.out_edges[target_node.position] = edge
            self.out_edges.append(edge)
            target_node.in_edges.append(edge)

        # Next_col must have sorted indices.
        assert all(n1.position < n2.position for n1, n2
                   in zip(next_col, next_col[1:]))

        # Iterate on column, calculating energy and adding edges as needed.
        col_iter = iter(next_col)

        # Find a node such that the index is greater than this node's.
        # TODO: account for maxBL and minBL
        curr_node = next(col_iter)
        while curr_node.position <= self.position:
            curr_node = next(col_iter)

        # Get the energy of the first valid node in the column, make edge.
        e = energy_func.block_energy(self.position, curr_node.position)
        add_edge(e, curr_node)

        # Start iterating over breakpoint columns starting at next position.
        start_position = curr_node.position + 1

        # Get the next node. If this is a StopIteration, we are at the end of
        # the column, so we can return.
        try:
            curr_node = next(col_iter)
        except StopIteration:
            return

        for i in range(start_position, next_col[-1].position):
            # Calculate energy difference using previous energy and energy
            # difference between previous and current indicies.
            e += energy_func.increment_block_energy(self.position, i)

            # If we're at curr_node's position, we can add an edge to the node.
            if i == curr_node.position:
                add_edge(e, curr_node)
                curr_node = next(col_iter)

        # Add the last edge. Doing it here avoids the StopIteration.
        add_edge(e, curr_node)

    def __hash__(self):
        return hash(repr(self))

    def __repr__(self):
        in_edges = [e.in_node.position for e in self.in_edges]
        out_edges = [e.out_node.position for e in self.out_edges]
        ret = f'Node({self.col}, {self.position}, {self.bp}, ' \
            f'in={in_edges}, out={out_edges})'
        return ret


class LibrariesNotFound(ValueError):
    pass


# TODO: ABC for optimizers like RASPP?


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

    Parameters:
        parents: Aligned parent sequences.
        n: Number of breakpoints in output libraries, not including breakpoints
            at positions 0 and len(parents.alignment).
        start_overhangs: Overhang options for the breakpoint at position 0. If
            None, no BreakPoint will be inserted at this location.
        start_overhangs: Overhang options for the breakpoint at position
            len(parents.alignment). If None, no BreakPoint will be inserted at
            this location.
        energy_func_type: Class used to calculate the energy of generated
            libraries.
        gg_enzyme: The restriction enzyme used to detect and evaluate Golden
            Gate sites.
        gg_threshold: The minimum threshold for a set of Golden Gate overhangs
            to be considered optimal. Setting this value lower may result in a
            faster runtime.
        amino_to_cdn: Mapping from amino acids to available codons. If None, a
            simple E. coli codon optimization will be used.

    Attributes:
        valid_bps (list[BreakPoint]): Parent alignment positions where a Golden
            Gate site could go. Libraries are constructed from combinations of
            elements in this list.
        columns (list[list[_Node]]): RASPP graph, e.g. figure 2 in Endelman et
            al. 2004. Traversal over this graph results in Libraries.

    """
    def __init__(
        self,
        parents: ParentSequences,
        n: int,
        start_overhangs: Optional[list[Overhang]] = None,
        end_overhangs: Optional[list[Overhang]] = None,
        energy_func_type: Type[EnergyFunction] = SCHEMA,
        gg_enzyme: Union[str, RestrictionEnzyme] = 'BsaI-HFv2',
        gg_threshold: Union[float, Decimal] = 0.95,
        amino_to_cdn: Optional[dict[str, set[str]]] = None,
    ):
        if isinstance(gg_enzyme, str):
            gg_enzyme = RestrictionEnzyme.from_name(gg_enzyme)

        energy_func = energy_func_type(parents)
        self._lib_config = LibraryConfig(energy_func, gg_enzyme, gg_threshold,
                                         None, amino_to_cdn)

        # Calculate all crossover sites valid for Golden Gate Assembly.
        self.valid_bps = calculate_breakpoints(
            parents,
            self._lib_config.amino_to_cdn,
            start_overhangs,
            end_overhangs
        )

        mr_cache = MutationRateCache.from_parents(parents, self.valid_bps)
        self._lib_config.mr_cache = mr_cache

        # Get the first column of RASPP graph.
        first_bp = sorted(self.valid_bps, key=attrgetter('position'))[0]
        if first_bp.position != 0:
            first_bp = BreakPoint(0, [])
        first_column = [_Node(0, first_bp)]

        # Build RASPP graph nodes. One column for each Library breakpoint plus
        # the starting column.
        columns = [first_column] + [[] for _ in range(n)]
        for pci, (prev_col, curr_col) in enumerate(zip(columns, columns[1:])):
            cci = pci + 1  # curr_col index from prev_col index
            min_pos = prev_col[0].position
            for bp in self.valid_bps:
                # Can exclude < min_pos because no edges will go into those.
                if bp.position >= min_pos:
                    node = _Node(cci, bp)
                    curr_col.append(node)

        # Fill out edges between nodes.
        columns[0][0].fill_out(columns[1], energy_func)
        for col, next_col in zip(columns[1:], columns[2:]):
            for node in col[:-1]:
                node.fill_out(next_col, energy_func)

        self.columns = columns

    def eval_bps(self, bps: Iterable[BreakPoint]) -> Decimal:
        """Average SCHEMA energy of sequence library with given breakpoints."""

        e = Decimal(0.0)
        curr_node = self.columns[0][0]
        for bp in bps:
            edge = curr_node.out_edges[bp.position]
            e += edge.e_diff
            curr_node = edge.out_node
        return e

    def _eval_all_bps(self) -> dict[tuple[BreakPoint, ...], Decimal]:
        """eval_bps over all combinations of breakpoints.

        Experimental method. TODO: Update or remove.

        """

        valid_bps = self.valid_bps.copy()
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
        min_blk_len: Optional[int] = None,
        max_blk_len: Optional[int] = None
    ) -> 'Library':
        """Find Library with minimum energy given block length constraints.

        Parameters:
            min_blk_len: Minimum block length. If None passed in, blocks can be
                arbitrarily small.
            max_blk_len: Maximum block length. If None passed in, blocks can be
                arbitrarily large.

        """

        n = len(self.columns) - 1  # number of breakpoints
        N = len(self._lib_config.energy_function.parents.alignment)

        if min_blk_len is None:
            min_blk_len = 0
        if max_blk_len is None:
            max_blk_len = N

        # Memoized traversal over columns of graph to find minimum library.
        # Record the minimum library for each node during traversal.
        first_node = self.columns[0][0]
        first_lib = _TempLibrary(Decimal(0.0), [first_node.bp])
        next_col_min_libs = {first_node: first_lib}
        while next_col_min_libs:
            curr_col_min_libs = next_col_min_libs

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
                    next_block_len = new_node.position - curr_node.position
                    if not min_blk_len <= next_block_len <= max_blk_len:
                        continue

                    # Get the best library for new_node so far and compare to
                    # library from curr_node. If no best library or curr_node
                    # library is best, set as the best library.
                    best_lib = next_col_min_libs.get(new_node)
                    if best_lib is None or curr_min_lib.energy + edge.e_diff \
                       < best_lib.energy:
                        new_lib = curr_min_lib.expand(edge.e_diff,
                                                      new_node.bp)
                        next_col_min_libs[new_node] = new_lib

        candidate_libs = [lib for node, lib in curr_col_min_libs.items()
                          if min_blk_len <= N - node.position <= max_blk_len]
        if not candidate_libs:
            raise LibrariesNotFound
        return_temp_lib = min(candidate_libs, key=lambda lib: lib.energy)
        if not return_temp_lib.breakpoints:
            raise LibrariesNotFound

        # Add last breakpoint, if necessary.
        last_bp = sorted(self.valid_bps)[-1]
        if last_bp.position == N:
            return_temp_lib = return_temp_lib.expand(Decimal(0.0), last_bp)

        return return_temp_lib.to_library(self._lib_config)

    def optimize(self, L_min: int, L_max: int) -> list[Library]:
        """Find energy-optimal libraries over a range of block lengths.

        Parameters:
            L_min: Minimum block length used for search.
            L_max: Maximum block length used for search.

        """

        chosen_libs = []

        # TODO: Can probably do problem reduction for each min_block_len
        # itethe cineration. E.g. For a given min_block_len iteration,
        # iterating
        # through max_block_len backwards allows us to say the most recently
        # found library is the minimum E until its max_block_len is greater
        # max_block_len.
        # TODO: DOUBLE CHECK THIS
        for min_block_len in range(L_min, L_max):
            curr_best = None
            for max_block_len in range(L_max, min_block_len, -1):
                # print(min_block_len, max_block_len, '\r', end='')
                if curr_best is not None:
                    if curr_best.max_block_len <= max_block_len:
                        continue
                    else:
                        curr_best = None
                try:
                    lib = self.min_bps(min_block_len, max_block_len)
                    curr_best = lib

                    if lib not in chosen_libs:
                        chosen_libs.append(lib)
                except LibrariesNotFound:
                    pass

        return chosen_libs

    def _all_libraries(
        self,
        min_block_len: Optional[int] = None,
        max_block_len: Optional[int] = None
    ) -> list['Library']:
        """Get all possible libraries with depth-first graph traversal.

        Experimental method. TODO: Update or remove.

        """
        libs = []
        n = len(self.columns) - 1  # number of breakpoints
        N = len(self._lib_config.energy_function.parents.alignment)

        if min_block_len is None:
            min_block_len = 0
        if max_block_len is None:
            max_block_len = N

        first_node = self.columns[0][0]
        first_lib = _TempLibrary(Decimal(0.0), [first_node.bp])
        stack = [(first_node, first_lib)]
        while stack:
            curr_node, curr_lib = stack.pop()
            if curr_node.col == n:
                next_block_len = N - curr_node.position
                if min_block_len <= next_block_len <= max_block_len:
                    last_bp = sorted(self.valid_bps,
                                     key=attrgetter('position'))[-1]
                    curr_temp_lib = curr_lib.expand(Decimal(0.0), last_bp)
                    e = curr_temp_lib.e
                    bps = curr_temp_lib.breakpoints
                    lib = Library.calc_from_config(bps, e, self._lib_config)
                    libs.append(lib)
                continue
            for edge in reversed(curr_node.out_edges):
                next_block_len = edge.out_node.position - curr_node.position
                if min_block_len <= next_block_len <= max_block_len:
                    new_node = edge.out_node
                    new_lib = curr_lib.expand(edge.e_diff, new_node.bp)
                    stack.append((new_node, new_lib))
        if not libs:
            msg = f'No libraries found with min_block_len={min_block_len}' \
                  f' and max_block_len={max_block_len}.'
            raise LibrariesNotFound(msg)
        return libs


def _generate_libraries(
    parents: ParentSequences,
    num_blocks: int,
    start_overhangs: Optional[list[Overhang]] = None,
    end_overhangs: Optional[list[Overhang]] = None,
    min_block_len: Optional[int] = None,
    max_block_len: Optional[int] = None,
    algorithm: str = 'SCHEMA-RASPP',
) -> list[Library]:
    """

    Wrapper for optimizers that generate chimeric protein libraries. Uses
    the restriction enzyme BsaI-HFv2, an initial Golden Gate efficiency
    threshold of 95%, and a simple E. coli codon optimization dictionary. For
    fine control over these parameters, directly use the Optimizer objects in
    the optimizer module.

    At present, the only implemented algorithm is SCHEMA-RASPP, which uses
    :class:`~schemarecomb.optimizers.RASPP` with the
    :class:`~schemarecomb.energy_functions.SCHEMA` energy function.

    The start_overhangs and end_overhangs parameters are commonly used to
    insert the assembled chimeras into a vector.

    Parameters:
        parents: Aligned parent sequences.
        num_blocks: Number of blocks in the generated libraries. Must be
            greater than 1.
        start_overhangs: Overhang options for the breakpoint at position 0. If
            None, no breakpoint will be inserted at this position.
        end_overhangs: Overhang options for the breakpoint at position
            len(parents.alignment). If None, no breakpoint will be inserted at
            this position.
        min_block_len: Smallest block length allowed in generated libraries. If
            None, the minimum block length will be len(parents.alignment) //
            (num_blocks + 1).
        max_block_len: Largest block length allowed in generated libraries. If
            None, the maximum block length will be len(parents.alignment) //
            (num_blocks - 1).
        algorithm: Name of the algorithm used to generate libraries. Currently
            must be 'SCHEMA-RASPP'.

    Returns:
        Collection of libraries found using the algorithm specified.

    Raises:
        NotImplementedError: If algorithm is anything except 'SCHEMA-RASPP'.
            This will be changed in future versions.
        ValueError: If num_blocks is less than 2.
        LibrariesNotFound: If no libraries could be generated given the inputs.

    """

    if algorithm != 'SCHEMA-RASPP':
        raise NotImplementedError('SCHEMA-RASPP is the only option available '
                                  'for the algorithm parameter at this time. ')
    if num_blocks <= 1:
        raise ValueError('num_blocks must be greater than 1.')

    # dict[str, tuple[type[EnergyFunction], type[Optimizer]]]
    algo_options = {'SCHEMA-RASPP': (SCHEMA, RASPP)}

    e_func_type, opt_type = algo_options[algorithm]

    opt = opt_type(
        parents,
        num_blocks-1,
        start_overhangs,
        end_overhangs,
        e_func_type
    )

    if min_block_len is None:
        min_block_len = len(parents.alignment) // (num_blocks + 1)
    if max_block_len is None:
        max_block_len = len(parents.alignment) // (num_blocks - 1)

    return opt.optimize(min_block_len, max_block_len)
