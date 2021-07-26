from dataclasses import dataclass, field
from functools import cached_property
from itertools import product, combinations

from ggsr.parent_alignment import ParentAlignment


@dataclass(repr=False)
class Library:
    energy: float
    breakpoints: dict[int, list[tuple[int, str]]] = field(default_factory=dict)

    # TODO: remove
    '''
    nonredun_bp_map: dict[int, int] = field(default_factory=dict)

    @property
    def nonredundant_breakpoint_indices(self):
        return tuple(self.nonredun_bp_map[bpi] for bpi in self.breakpoints)
    '''

    @cached_property
    def min_block_len(self):
        bps = sorted(self.breakpoints)
        return min(bp2 - bp1 for bp1, bp2 in zip(bps, bps[1:]))

    @cached_property
    def max_block_len(self):
        bps = sorted(self.breakpoints)
        return max(bp2 - bp1 for bp1, bp2 in zip(bps, bps[1:]))

    def expand(self, e_diff: float, new_bp: dict[int, list[tuple[int, str]]]):
        new_e = self.energy + e_diff
        new_bps = self.breakpoints | new_bp
        return type(self)(new_e, new_bps)

    def __repr__(self):
        return f'Library({self.energy}, {list(self.breakpoints.keys())})'

    @property
    def average_m(self):
        return self._m

    def calc_average_m_naive(self, pa: ParentAlignment, bp_to_group,
                             group_bp_cache):
        """Calculate the average number of mutations in library.

        Naive way to calculate <M>. Currently kept for testing purposes, may be
        deleted later for simplicity.
        """
        def calc_muts(seq1, seq2):
            return sum(1 for a1, a2 in zip(seq1, seq2) if a1 != a2)

        bps = tuple(self.breakpoints)

        # ADDITION
        bp_groups = tuple(bp_to_group[bp] for bp in bps)
        if (m := group_bp_cache.get(bp_groups)) is not None:
            self._m = m
            return m
        # END ADDITION

        p_seqs = [str(sr.seq) for sr in pa.aligned_sequences]

        p = len(p_seqs)  # number of parents
        N = len(p_seqs[0])  # number of amino acids in alignment

        # bps must be (0, ..., N).
        if bps[0] != 0:
            bps = (0,) + bps
        if bps[-1] != N:
            bps += (N,)

        n = len(bps) - 2  # number of crossovers. number of blocks is n+1

        # Construct all library chimeras to calculate average mutations.
        total_muts = 0
        for i, block_parents in enumerate(product(p_seqs, repeat=n+1)):
            # Construct chimera from parent blocks.
            block_seqs = [blkpar[start: end] for blkpar, start, end
                          in zip(block_parents, bps, bps[1:])]
            chimera = ''.join(block_seqs)

            # Chimera m is the number of mutations from the closest parent.
            total_muts += min(calc_muts(chimera, parent) for parent in p_seqs)

        m = total_muts / p**(n+1)

        # ADDITION
        # global group_bp_map
        group_bp_cache[bp_groups] = m
        # END ADDITION

        self._m = m
        return m

    def calc_average_m(self, pa: ParentAlignment, bp_to_group, group_bp_cache):
        """Calculate the average number of mutations in library."""
        def sequence_mutations(seq1, seq2):
            """Hamming distance between seq1 and seq2."""
            assert len(seq1) == len(seq2)
            return sum(1 for a1, a2 in zip(seq1, seq2) if a1 != a2)

        bp_groups = tuple(bp_to_group[bp] for bp in self.breakpoints)
        if (M := group_bp_cache.get(bp_groups)) is not None:
            self._m = M
            return M

        seqs = [str(sr.seq) for sr in pa.aligned_sequences]
        blmuts = []
        bps = sorted(self.breakpoints)
        for bp1, bp2 in zip(bps, bps[1:]):
            muts = {i: {i: 0} for i, _ in enumerate(seqs)}
            for (i, s1), (j, s2) in combinations(enumerate(seqs), 2):
                s1_slice = s1[bp1:bp2]
                s2_slice = s2[bp1:bp2]
                m = sequence_mutations(s1_slice, s2_slice)
                muts[i][j] = m
                muts[j][i] = m
            blmuts.append(muts)

        # Calculate parameters for mutation calculation.
        num_blocks = len(self.breakpoints) - 1
        num_parents = len(seqs)

        # Perform depth first traversal over blocks to calculate M_sum. This
        # is faster than itertools.product because it avoids redundant
        # summation.
        stack = [(1, [blmuts[0][p1][p2] for p2 in range(num_parents)])
                 for p1 in range(num_parents)]
        M_sum = 0
        while stack:
            i, par_muts = stack.pop()
            par_iter = list(zip(par_muts, range(num_parents)))
            block_muts = blmuts[i]
            for p1 in range(num_parents):
                p_muts = block_muts[p1]
                new_part_muts = [m + p_muts[p2] for m, p2 in par_iter]
                if i + 1 == num_blocks:
                    M_sum += min(new_part_muts)
                else:
                    stack.append((i + 1, new_part_muts))

        M = M_sum / num_parents**num_blocks
        self._m = M

        group_bp_cache[bp_groups] = M

        return M

    '''
    def calc_average_m_rust(self, pa: ParentAlignment, bp_to_group):
        import ggsr_rust

        def sequence_mutations(seq1, seq2):
            return sum(1 for a1, a2 in zip(seq1, seq2) if a1 != a2)

        seqs = [str(sr.seq) for sr in pa.aligned_sequences]
        blmuts = []
        bps = sorted(self.breakpoints)
        for bp1, bp2 in zip(bps, bps[1:]):
            # muts = {i: {i: 0} for i, _ in enumerate(seqs)}
            muts = np.empty((len(seqs), len(seqs)), dtype=int)
            np.fill_diagonal(muts, 0)
            for (i, s1), (j, s2) in combinations(enumerate(seqs), 2):
                s1_slice = s1[bp1:bp2]
                s2_slice = s2[bp1:bp2]
                m = sequence_mutations(s1_slice, s2_slice)
                muts[i, j] = m
                muts[j, i] = m
            blmuts.append(muts.tolist())

        M = ggsr_rust.calc_average_m(blmuts)
        self._m = M

        return M
    '''
