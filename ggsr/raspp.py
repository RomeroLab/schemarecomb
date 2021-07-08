"""Module for running RASPP algorithm."""

from collections import namedtuple
from dataclasses import dataclass
from functools import cached_property
from itertools import product, combinations

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
    def __init__(self, pa):
        alignment = list(zip(*pa.aligned_sequences))
        self.E_matrix = np.zeros((len(alignment), len(alignment)))
        for i, j in pa.pdb_structure.contacts:
            AAs_i, AAs_j = alignment[i], alignment[j]
            parental = set(zip(AAs_i, AAs_j))
            broken = sum(1 for combo in product(AAs_i, AAs_j) if combo
                         not in parental)
            self.E_matrix[i, j] = broken / (len(AAs_i)**2)

    def block_energy(self, start, end):
        """ Calculate the SCHEMA energy of a block.

        Args:
            start: Int block's start breakpoint index.
            end: Int block's end breakpoint index.
            E_matrix: SCHEMA energy matrix from generate_weighted_E_matrix.
        """
        energy = self.E_matrix[start:end, end:].sum()
        return energy

    def diff(self, x_k_1, x_k_inc):
        x_k = x_k_inc - 1  # previous x_k
        neg = self.E_matrix[x_k_1:x_k, x_k].sum()
        pos = self.E_matrix[x_k, x_k_inc:].sum()
        return pos - neg


class RASPP:
    def __init__(self, parent_aln: ParentAlignment,
                 run_params: dict,  # change to class later
                 cdn_candidates: dict = AA_C31,  # change later to class?
                 energy_f: type[EnergyFunction] = SCHEMA):
        self.aln = parent_aln
        self.cdns = cdn_candidates
        self.e = SCHEMA(parent_aln)
        self.vector_overhangs = run_params['vector_overhangs']

        self.columns = [[Node(0, 0, [], {})]]
        N = len(parent_aln.aligned_sequences[0])
        for col in range(1, run_params['n'] + 1):
            if run_params['all_nodes']:
                # prev col will never use n.index less than col
                new_col = [Node(col, index, [], {}) for index in range(col, N)]
            else:
                # prev col will never use n.index less than col
                col_bps = [index for index in self.breakpoints if index >= col]
                new_col = [Node(col, index, [], {}) for index in col_bps]
            self.columns.append(new_col)

        for col, n_col in zip(self.columns, self.columns[1:]):
            if len(col) == 1:
                col[0].fill_out(self.aln, n_col, self.e.block_energy,
                                self.e.diff)
            else:
                for node in col[:-1]:
                    node.fill_out(self.aln, n_col, self.e.block_energy,
                                  self.e.diff)

    @cached_property
    def breakpoints(self):
        vector_overhangs = self.vector_overhangs
        aln_seqs = [str(sr.seq) for sr in self.aln.aligned_sequences]

        # Amino acids at each position.
        aln_sets = zip(*aln_seqs)

        # Sets of codons for each amino acid at an alignment position.
        aln_cdns = [{self.cdns[aa] for aa in pos_aas} for pos_aas in aln_sets]

        breakpoints = {0: [vector_overhangs[0]]}

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

        breakpoints[len(aln_cdns)] = [vector_overhangs[1]]

        # TODO: Check if there's enough breakpoints for the run params.

        return breakpoints

    def eval_bps(self, bps):
        e = 0.0
        curr_node = self.columns[0][0]
        for bp in bps:
            edge = curr_node.edges[bp]
            e += edge.e_diff
            curr_node = edge.out_node
        return e

    def eval_all_bps(self):
        valid_bps = list(self.breakpoints.keys())
        try:
            del valid_bps[0]
            del valid_bps[-1]
        except KeyError:
            pass
        for bps in combinations(valid_bps, len(self.columns)-1):
            self.eval_bps(bps)

    def max_bps(self):
        col_iter = iter(self.columns)

        curr_col = next(col_iter)
        curr_col[0].min_eng = 0.0
        curr_col[0].min_bps = tuple([])

        for col in col_iter:
            for node in col:
                node.min_eng = 100000000000  # large number
                if not node.in_edges:
                    continue
                for edge in node.in_edges:
                    if not hasattr(edge.in_node, 'min_bps'):
                        continue
                    new_eng = edge.in_node.min_eng + edge.e_diff
                    if new_eng < node.min_eng:
                        node.min_eng = new_eng
                        node.min_bps = edge.in_node.min_bps + (node.index,)
                if not hasattr(node, 'min_bps'):
                    continue


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


# might want to add edge to this
Edge = namedtuple('Edge', ['e_diff', 'in_node', 'out_node'])


@dataclass(repr=False)
class Node:
    col: int
    index: int
    in_edges: list
    edges: dict

    def fill_out(self, pa, next_col, eng, eng_diff):
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
        curr_node = next(col_iter)
        while curr_node.index <= self.index:
            curr_node = next(col_iter)

        # Get the energy of the first valid node in the column, make edge.
        e = eng(self.index, curr_node.index)
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
            e += eng_diff(self.index, i)

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
