'''Library for choosing additional sequences for parent alignments.

choose_candidates is the main function. The TreeNode and Tree classes help
choose_candidates implement a tree-based search over all possible length-
<num_additional> combinations of candidate sequences. Importantly, the search
is short circuited such that the function returns if the best combination is
known to be found early.
'''

from .utils import _calc_identity


class TreeNode:
    def __init__(self, cands=[], max_q=0.0, max_pq=0.0, parent=None):
        self.cands = cands  # list of candidate SRs
        self.children = []  # Tree children, ordered by max_q
        self.max_q = max_q  # max diff between cands
        self.max_pq = max_pq  # max diff between cands and existing parents
        self.parent = parent

    def new_cand(self, cand, cand_diff, diff_thresh):
        '''Create new node with self.cands + [cand]'''
        if self.cands:
            new_q = max(abs(0.7 - _calc_identity(cand, x)) for x in self.cands)
            new_q = max(self.max_q, new_q)
        else:
            new_q = 0.0
        if max(new_q, cand_diff) >= diff_thresh:
            # don't create a new node because it is already non-optimal
            return None  # probably change b/c EP2 Item 20
        new_cands = self.cands + [cand]
        return TreeNode(new_cands, new_q, cand_diff, self)

    def add_cand(self, cand, cand_diff, diff_thresh):
        '''Create new child node with self.cands + [cand]'''
        new_node = self.new_cand(cand, cand_diff, diff_thresh)
        if new_node is not None:
            self.children.append(new_node)

    def delete(self):
        self.parent.children.remove(self)

    def __repr__(self):
        return ','.join([c.id for c in self.cands])


class Tree:
    def __init__(self, num_seqs):
        self.base = TreeNode()
        self.num_seqs = num_seqs  # number of sequences to choose
        self.best_leaf = None
        self.best_diff = 1.0

    def add_cand(self, cand, cand_diff):
        # recursive breadth-first
        stack = [self.base]
        while stack:
            curr_node = stack.pop()

            # prune
            curr_max = max(curr_node.max_q, curr_node.max_pq)
            if curr_max > self.best_diff:
                curr_node.delete()
                continue

            # check if leaf
            if len(curr_node.cands) == self.num_seqs - 1:
                leaf = curr_node.new_cand(cand, cand_diff, self.best_diff)
                if leaf is not None:
                    # leaf max diff is smaller than self.best_diff
                    self.best_leaf = leaf
                    leaf_diff = max(leaf.max_q, leaf.max_pq)
                    self.best_diff = leaf_diff
                    if leaf.max_pq >= leaf.max_q:
                        return 'best found'
                continue

            stack += curr_node.children  # comes first b/c don't want dups
            curr_node.add_cand(cand, cand_diff, self.best_diff)

    def shape(self):
        from collections import defaultdict
        shape = defaultdict(int)
        stack = [self.base]
        while stack:
            curr_node = stack.pop()
            shape[len(curr_node.cands)] += 1
            stack += curr_node.children
        return [shape[i] for i, _ in enumerate(shape)]


def choose_candidates(candidate_sequences, existing_parents, num_additional,
                      desired_identity):

    # For each candidate, find the maximum identity difference with each
    # parent.
    cand_diffs = []
    for i, cand in enumerate(candidate_sequences):
        print(i, '\r', end='')
        max_diff = max(abs(desired_identity - _calc_identity(cand, x))
                       for x in existing_parents)
        cand_diffs.append((cand, max_diff))

    # Sort candidates by difference to enable short-circuiting.
    # sorted is fast enough (but could be faster with heapq)
    sorted_cand_diffs = list(sorted(cand_diffs, key=lambda x: x[1]))

    # TODO: assure sorted_cand_diffs is sorted
    tree = Tree(num_additional)
    for i, (cand, cand_diff) in enumerate(sorted_cand_diffs):
        ret = tree.add_cand(cand, cand_diff)
        if ret == 'best found':
            break
        if i > 5:
            break
    return tree.best_leaf.cands
