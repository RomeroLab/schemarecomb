from Bio import pairwise2


def _calc_identity(seq1, seq2):
    aln_result, = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    return aln_result.score / (aln_result.end + 1)


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
        print(self.best_diff)
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


def minmax(sorted_cand_diffs, num_additional):
    # TODO: assure sorted_cand_diffs is sorted
    print(num_additional)
    tree = Tree(num_additional)
    print(tree.shape())
    for i, (cand, cand_diff) in enumerate(sorted_cand_diffs):
        ret = tree.add_cand(cand, cand_diff)
        print(i, tree.shape())
        if ret == 'best found':
            break
        if i > 5:
            break
    print(tree.best_diff)
    print(tree.best_leaf)
    return tree.best_leaf.cands
