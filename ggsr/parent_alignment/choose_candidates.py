# choose_candidates.py

"""Library for choosing additional sequences for parent alignments.

choose_candidates is the main function. The TreeNode and Tree classes help
choose_candidates implement a tree-based search over all possible length-
<num_additional> combinations of candidate sequences. Importantly, the search
is short circuited such that the function returns if the best combination is
known to be found early.
"""

from collections import defaultdict

from Bio import SeqRecord

from .utils import _calc_identity


class TreeNode:
    """Represents a set of candidate sequences.

    Attributes:
        cands: List of candidates SeqRecords.
        max_cc_diff: Maximum out of |% identity - target_identity| between
        each pair of candidates.
        max_pc_diff: Maximum out of |% identity - target_identity| of each
            parent and each candidate.
        children: Child nodes of this node.
        parent: Parent node of this node.
    """
    def __init__(self, cands: SeqRecord.SeqRecord = [],
                 max_cc_diff: float = 0.0,
                 max_pc_diff: float = 0.0,
                 parent: 'TreeNode' = None) -> None:
        """Initialize node instance with precomputed identies."""
        self.cands = cands
        self.max_cc_diff = max_cc_diff
        self.max_pc_diff = max_pc_diff
        self.children = []
        self.parent = parent

    def new_cand(self, cand: SeqRecord.SeqRecord, new_pc_diff: float,
                 diff_thresh: float, target_identity: float) -> 'TreeNode':
        """Create new child node with candidates self.cands + [cand].

        Args:
            cand: New candidate sequence.
            new_pc_diff: Maximum out of |% identity - target_identity| between
                cand and each parent. This is assumed to be greater than the
                "pc_diff" of anything in self.cands, by virtue of pre-sorting
                based on the pc_diff.
            diff_thresh: Current known optimum for a valid set of candidates.
                If the potentially new set has a higher diff than this value,
                we can skip creation because this node and all its children
                are non-optimal.
            target_identity: Ideal cross-wise identity between all sequences in
                the concatenation set of parents and selected candidates.
        """
        if self.cands:
            # Calculate the identity diff between the new candidate and each
            # candidate in self.cands.
            new_cc_diff = max(abs(target_identity - _calc_identity(cand, x))
                              for x in self.cands)
            # Choose maximum between old cc_diff and new_cc_diff.
            new_cc_diff = max(self.max_cc_diff, new_cc_diff)
        else:
            new_cc_diff = 0.0  # only one candidate in new node
        if max(new_cc_diff, new_pc_diff) >= diff_thresh:
            # don't create a new node because it is already non-optimal
            return None  # probably change b/c EP2 Item 20
        new_cands = self.cands + [cand]
        return TreeNode(new_cands, new_cc_diff, new_pc_diff, self)

    def add_cand(self, cand: SeqRecord.SeqRecord, new_pc_diff: float,
                 diff_thresh: float, target_identity: float) -> None:
        """new_cand new_cand for adding non-leaf nodes."""
        new_node = self.new_cand(cand, new_pc_diff, diff_thresh,
                                 target_identity)
        if new_node is not None:
            self.children.append(new_node)

    def delete(self) -> None:
        """Delete node from tree."""
        self.parent.children.remove(self)

    def __repr__(self) -> str:
        """Return CSV string of the node's candidate ids."""
        return ','.join([c.id for c in self.cands])


class Tree:
    """Data structure for finding the ideal set of candidate sequences.

    Attributes:
        base: Node representing the empty set of candidates.
        num_seq: Number of sequences to choose. Depth of tree is num_seq + 1.
        target_identity: Ideal cross-wise identity between all sequences in
            the concatenation set of parents and selected candidates.
        best_leaf: The TreeNode with the current best set of candidates.
        best_diff: The maximum |% identity - target_identity| out of each pair
            in best_leaf.cands + parents.
    """
    def __init__(self, num_seqs: int, target_identity: float) -> None:
        """Tree for finding <num_seqs> candidates with <target_identity>."""
        self.base = TreeNode()
        self.num_seqs = num_seqs
        self.target_identity = target_identity
        self.best_leaf = None
        self.best_diff = 1.0

    def add_cand(self, cand: SeqRecord.SeqRecord, new_pc_diff: float) -> None:
        """Add cand to candidate sets in tree and add new nodes if needed.

        Args:
            cand: New candidate sequence.
            new_pc_diff: Maximum out of |% identity - target_identity| between
                cand and each parent. This is assumed to be greater than the
                "pc_diff" of anything in self.cands, by virtue of pre-sorting
                based on the pc_diff.
        """
        # Breadth-first traversal over tree, adding new cand to each node.
        stack = [self.base]
        while stack:
            curr_node = stack.pop()

            # Prune tree if curr_node's max identity difference is greater
            # than the current known best.
            curr_max = max(curr_node.max_cc_diff, curr_node.max_pc_diff)
            if curr_max > self.best_diff:
                curr_node.delete()
                continue

            # Leaf node case: if new node will have num_seqs cands, then we can
            # evaluate it directly and don't need to add it to tree.
            if len(curr_node.cands) == self.num_seqs - 1:
                leaf = curr_node.new_cand(cand, new_pc_diff, self.best_diff,
                                          self.target_identity)
                if leaf is not None:
                    # leaf max diff is smaller than self.best_diff, new best!
                    self.best_leaf = leaf
                    self.best_diff = max(leaf.max_cc_diff, leaf.max_pc_diff)

                    # Every future leaf will have equal or greater max_pc_diff,
                    # so this leaf must be the best one and we can stop.
                    if leaf.max_pc_diff >= leaf.max_cc_diff:
                        return 'best found'

                continue

            # Non-leaf node case.
            stack += curr_node.children  # comes first b/c don't want dup cands
            curr_node.add_cand(cand, new_pc_diff, self.best_diff,
                               self.target_identity)

    def shape(self) -> list[int]:
        """Number of nodes at each level in tree."""
        shape = defaultdict(int)
        stack = [self.base]
        while stack:
            curr_node = stack.pop()
            shape[len(curr_node.cands)] += 1
            stack += curr_node.children
        return [shape[i] for i, _ in enumerate(shape)]


def choose_candidates(candidate_sequences: list[SeqRecord.SeqRecord],
                      existing_parents: list[SeqRecord.SeqRecord],
                      num_additional: int,
                      desired_identity: float) -> list[SeqRecord.SeqRecord]:
    """Choose the ideal set of candidate sequences.

    Ideal set is defined as the set of candidates with the minimum max_diff,
    where max_diff is the maximum cross-wise |% identity - desired_identity|
    between each pair of sequences in the concatenation of the set and parents.

    Args:
        candidate_sequences: Sequences able to be selected.
        existing_parents: Parent sequences in the library already.
        num_additional: Number of candidate sequences to choose.
        desired_identity: Ideal cross-wise identity between all sequences in
            the concatenation set of parents and selected candidates.

    Returns:
        Ideal set of candidate sequences.
    """

    # For each candidate, find the maximum identity difference with each
    # parent.

    cand_diffs = []
    for i, cand in enumerate(candidate_sequences):
        print(i, '\r', end='')
        max_diff = max(abs(desired_identity - _calc_identity(cand, x))
                       for x in existing_parents)
        cand_diffs.append((cand, max_diff))
    print()

    # Sort candidates by difference to enable short-circuiting. sorted is fast
    # enough, but this could be faster with heapq.
    sorted_cand_diffs = list(sorted(cand_diffs, key=lambda x: x[1]))

    # Construct Tree and find the best set. In the worst case this might take
    # awhile.
    # TODO: Bound the time this requires?
    tree = Tree(num_additional, desired_identity)
    for cand, pc_diff in sorted_cand_diffs:
        ret = tree.add_cand(cand, pc_diff)
        if ret == 'best found':
            break
    return tree.best_leaf.cands
