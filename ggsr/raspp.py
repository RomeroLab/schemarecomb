"""Module for running the RASPP algorithm.

TODO: Refactor to separate codon_options, EnergyFunctions.
"""

from collections.abc import Iterable
from collections import namedtuple
from dataclasses import dataclass
from itertools import product, combinations
from typing import Optional, Union

import numpy as np

from ggsr.parent_alignment import ParentAlignment
from ggsr.parent_alignment.pdb_structure import PDBStructure


# List of optimal codons for E. coli.
AA_C31 = {'A': ('GCT', 'GCA'), 'R': ('CGT', 'CGA'), 'N': ('AAT',),
          'D': ('GAT',), 'C': ('TGT',), 'Q': ('CAA', 'CAG'), 'E': ('GAA',),
          'G': ('GGT',), 'H': ('CAT', 'CAC'), 'I': ('ATT', 'ATC'),
          'L': ('TTA', 'TTG', 'CTA'), 'K': ('AAA',), 'M': ('ATG',),
          'F': ('TTT',), 'P': ('CCT', 'CCA'), 'S': ('AGT', 'TCA'),
          'T': ('ACA', 'ACT'), 'W': ('TGG',), 'Y': ('TAT',),
          'V': ('GTT', 'GTA'), '*': ('TGA',), '-': ('---',)}


class EnergyFunction:
    pass


class SCHEMA(EnergyFunction):
    """SCHEMA energy.

    Given a PDB structure and a multiple sequence alignment, the SCHEMA energy
    of a recombinant sequence is the number of pairwise interactions broken by
    recombination. Expressed as an equation (Endelman et al. 2004 eq. 7):

        E = sum_i sum_{j>i} C_{i,j} D_{i,j}

        where:
            C_{i,j} = 1 if residues i and j are within 4.5 Angstroms, else 0.
            D_{i,j} = 0 if the amino acids at positions i and j in the
                recombinant sequence are also present together in any parent,
                else 0.

    Example: If "AEH" and "GQH" are aligned parents with contacts (0, 1) and
        (1, 2), chimera "AQH" has a SCHEMA energy of 1.0 because residues 'A'
        and 'Q' are not found in any parents at positions 0 and 1,
        respectively, and there is a contact at these positions. Although
        residues 1 and 2 are contacting, 'Q' and 'H' are found together at
        these positions in the second parent, so D_{0,1} = 0 and doesn't
        contribute to the SCHEMA energy.

    See Voigt et al. 2002 and Silberg et al. 2004 for more information.

    Public methods:
        block_energy: Average SCHEMA energy of a section of a parent multiple
            sequence alignment.
        increment_block_energy: Average SCHEMA energy difference between a
            block and a block decreased by one residue.
    """
    def __init__(self, pa: ParentAlignment):
        """Generate the SCHEMA energy matrix for block calculation.

        E_matrix is a square matrix with the same size as the length of
        alignment. E_matrix_{r,j} is Endelman et al. 2004 eq. 6. with the r and
        t sums removed and evaluated for average SCHEMA energy for a specific r
        and t. Then eq. 6 can be evaluated for a given block by summing over
        E_matrix_{r,t}.

        Args:
            pa: Parent sequence alignment. Must have pdb_structure attribute.
        """
        alignment = list(zip(*pa.aligned_sequences))
        self.E_matrix = np.zeros((len(alignment), len(alignment)))

        # Only need to evaluate at contacts because otherwise C_{i,j} = 0.
        for i, j in pa.pdb_structure.contacts:
            # Get the ordered parent amino acids at i and j.
            AAs_i, AAs_j = alignment[i], alignment[j]

            # Parent amino acid pair might be redundant.
            parental = set(zip(AAs_i, AAs_j))

            # For each possible recombinant calculate D_ij. E_matrix entry
            # is this sum over the number of possible recombinants.
            broken = sum(1 for combo in product(AAs_i, AAs_j) if combo
                         not in parental)
            self.E_matrix[i, j] = broken / (len(AAs_i)**2)

    def block_energy(self, start: int, end: int) -> float:
        """Calculate the average SCHEMA energy of a block.

        This is Endelman et al. 2004 eq 6.

        Args:
            start: Index of the first residue in the block.
            end: Index of the first residue in the next block. Last residue in
                the current block is at end-1.

        Returns:
            Average SCHEMA energy of block [start, end].
        """
        energy = self.E_matrix[start:end, end:].sum()
        return energy

    def increment_block_energy(self, start: int, new_end: int) -> float:
        """Difference in average SCHEMA energy when block size is incremented.

        This is equivalent to
            block_energy(start, new_end) - block_energy(start, new_end-1).

        This function is faster than separate block_energy calculations.
        Therefore, you can calculate the energy for the first node in the
        graph using block_energy, then calculate the subsequent node energies
        using this function.

        Args:
            start: Index of the first residue in the block.
            end: Index of the first residue in the next block, after the
                current block size is incremented.

        Returns:
            Difference in average SCHEMA energy between block [start, new_end]
                and block [start, new_end-1].
        """
        old_end = new_end - 1
        neg = self.E_matrix[start:old_end, old_end].sum()
        pos = self.E_matrix[old_end, new_end:].sum()
        return pos - neg


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
            codon_options: dict[str: list[str]] = AA_C31,
            energy_func: Union[EnergyFunction, type[EnergyFunction]] = SCHEMA
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
        if isinstance(energy_func, type):
            energy_func = energy_func(parent_aln)

        # Calculate all crossover sites valid for Golden Gate Assembly.
        self.breakpoints = _calculate_breakpoints(parent_aln, codon_options,
                                                  start_overhang, end_overhang)

        # Build RASPP graph nodes.
        self.columns = [[Node(0, 0, [], {})]]
        for col in range(1, n + 1):
            # No edges into node with index less than col => exclude.
            col_indices = [index for index in self.breakpoints if index >= col]
            new_col = [Node(col, index, [], {}) for index in col_indices]
            self.columns.append(new_col)

        # Fill out edges between nodes.
        self.columns[0][0].fill_out(parent_aln, self.columns[1], energy_func)
        for col, next_col in zip(self.columns[1:], self.columns[2:]):
            for node in col[:-1]:
                node.fill_out(parent_aln, next_col, energy_func)

    def eval_bps(self, bps: Iterable[int]) -> float:
        """Average SCHEMA energy of sequence library with given breakpoints."""
        e = 0.0
        curr_node = self.columns[0][0]
        for bp in bps:
            edge = curr_node.edges[bp]
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

    def min_bps(self) -> tuple[tuple[int], float]:
        """Find the set of breakpoints with minimum energy."""
        col_iter = iter(self.columns)

        curr_col = next(col_iter)
        curr_col[0].min_energy = 0.0
        curr_col[0].min_bps = ()

        for col in col_iter:
            for node in col:
                node.min_energy = 100000000000  # large number
                for edge in node.in_edges:
                    new_energy = edge.in_node.min_energy + edge.e_diff
                    if new_energy < node.min_energy:
                        node.min_energy = new_energy
                        node.min_bps = edge.in_node.min_bps + (node.index,)

        best_node = min(self.columns[-1], key=lambda x: x.min_energy)
        return best_node.min_bps, best_node.min_energy


def _get_valid_patterns(cdn_sets, reverse=False):
    fwd_lim, bwd_lim = int(reverse), 1 - int(reverse)

    def cdn_reduce(cdns):
        return {cdn[fwd_lim:len(cdn)-bwd_lim] for cdn in cdns}

    # 1 AA per codon, so len 3 patterns exist <=> len(cdn_sets) == 1.
    if len(cdn_sets) == 1:
        len_3_sets = tuple(cdn_sets)[0]
    else:
        len_3_sets = {}

    cdn_sets = [cdn_reduce(cdns) for cdns in cdn_sets]
    len_2_sets = set.intersection(*cdn_sets)

    cdn_sets = [cdn_reduce(cdns) for cdns in cdn_sets]
    len_1_sets = set.intersection(*cdn_sets)

    return len_1_sets, len_2_sets, len_3_sets


def _calculate_breakpoints(parent_aln, codon_options, start_overhang,
                           end_overhang):
    aln_seqs = [str(sr.seq) for sr in parent_aln.aligned_sequences]

    # Amino acids at each position.
    aln_sets = zip(*aln_seqs)

    # Sets of codons for each amino acid at an alignment position.
    aln_cdns = [{codon_options[aa] for aa in pos_aas} for pos_aas in aln_sets]

    if start_overhang is not None:
        breakpoints = {0: [start_overhang]}

    # Search possible Golden Gate sites by iterating over adjacent codons.
    for bp, (cdns1, cdns2) in enumerate(zip(aln_cdns, aln_cdns[1:]), 1):
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
            continue

        breakpoints[bp] = overhangs

    if end_overhang is not None:
        breakpoints[len(aln_cdns)] = [end_overhang]

    # TODO: Check if there's enough breakpoints for the run params.

    return breakpoints


# might want to add edge to this
Edge = namedtuple('Edge', ['e_diff', 'in_node', 'out_node'])


@dataclass(repr=False)
class Node:
    col: int
    index: int
    in_edges: list
    edges: dict

    def fill_out(self, pa, next_col, energy_func):
        """Fill out edges from this node to next column."""

        def add_edge(energy, target_node):
            """Add an edge from self to target_node."""
            edge = Edge(e, self, target_node)
            self.edges[target_node.index] = edge
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


if __name__ == '__main__':
    '''
    loc = '../tests/bgl3_sample/truncated/'
    pa = ParentAlignment.from_fasta(loc+'trunc.fasta')
    pdb = PDBStructure.from_pdb_file(loc+'trunc.pdb')
    '''
    import time
    s = time.time()
    loc = '../tests/bgl3_sample/'
    pa = ParentAlignment.from_fasta(loc+'bgl3_sequences.fasta')
    pdb = PDBStructure.from_pdb_file(loc+'1GNX.pdb')
    pa.pdb_structure = pdb

    vector_overhangs = [(0, 'TATG'), (3, 'TGAG')]
    run_params = {'n': 2, 'all_nodes': False,
                  'vector_overhangs': vector_overhangs}
    # TODO: add minBL, maxBL
    s1 = time.time()
    raspp = RASPP(pa, run_params)
    r_time = time.time() - s1
    # raspp.eval_all_bps()
    s1 = time.time()
    raspp.max_bps()
    print('raspp', r_time)
    print('max', time.time() - s1)
    print(time.time() - s)
