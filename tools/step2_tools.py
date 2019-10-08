""" This module contains tools used by the Step2 script. """

import collections
import itertools
import os
import pickle
import re
import sys

import networkx
import numpy

from tools import general_tools

# Change these as needed
dirname = os.path.split(os.path.abspath(os.path.realpath(sys.argv[0])))[0]
lig_counts_fn = dirname + '/tools/normalized_ligation_counts_18h_37C.p'
normalized_ligation_counts = pickle.load(open(lig_counts_fn, 'rb'))
threshold = 245  # ligation count threshold for eliminating bad overhangs

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def complementary_sequence(sequence):
    """ Returns the DNA complement of sequence. """
    return ''.join(reversed([complement[bp] for bp in sequence]))


def acceptable_overhang(overhang, vector_overhangs):
    """Returns true if overhang is not palindromic, full-GC, below threshold,
    or is not a vector_overhang. """

    # palindromic overhangs
    if overhang == complementary_sequence(overhang):
        return False

    # high GC overhangs
    if 'A' not in overhang and 'T' not in overhang:
        return False

    # overhangs with low rates of correct ligations
    if normalized_ligation_counts[overhang][complementary_sequence(overhang)] \
       < threshold:
        return False

    # vector overhangs
    if overhang in vector_overhangs or complementary_sequence(overhang) \
       in vector_overhangs:
        return False

    return True


all_overhangs = [''.join(oh) for oh in itertools.product(general_tools.bases,
                                                         repeat=4)]


def _pattern_valid(pattern, pos_cdns):
    """ Check if pattern can represent pos_cdns. Helper function for
    get_valid_patterns.

    Args:
        pattern: str of length 3 that represents a codon pattern. Wild card
            base positions are denoted with '.' for re compatibility.
        pos_cdns: List of lists of codons. Each inner list has all valid codons
            for one of the unique residues at the given MSA position.

    Returns:
        True if at least one codon in each pos_cdns element follows pattern.
    """
    for AA_cdns in pos_cdns:
        matches_a_cdn = any(re.fullmatch(pattern, codon) for codon in AA_cdns)
        if not matches_a_cdn:
            return False
    return True


def _get_valid_patterns(pos_cdns, reverse=False):
    """ Gets patterns that describe at least one codon in each pos_cdn element.

    Each pattern is a length 3 str that may use '.' as a wild card.
    Patterns are explored left to right e.g. if 'A..' works then 'AA.', 'AG.'
    are explored, unless reverse=True.

    Args:
        pos_cdns: List of lists of codons. Each inner list has all valid codons
            for one of the unique residues at the given MSA position.
        reverse: Boolean specifying if patterns should be explored right to
            left.

    Returns:
        List of valid patterns.
    """
    # reverse cdns if pattern search is at end of cdns rather than start
    if reverse:
        pos_cdns = [[cdn[::-1] for cdn in AA_cdns] for AA_cdns in pos_cdns]

    # need to iterate through all possible cdn patterns, but can eliminate
    # based on subpattern validity
    valid_patterns = []
    for b1 in general_tools.bases:
        pattern = f'{b1}..'
        if not _pattern_valid(pattern, pos_cdns):
            continue  # No use exploring derivatives
        valid_patterns.append(pattern)
        for b2 in general_tools.bases:
            pattern = f'{b1}{b2}.'
            if not _pattern_valid(pattern, pos_cdns):
                continue  # No use exploring derivatives
            valid_patterns.append(pattern)
            for b3 in general_tools.bases:
                pattern = f'{b1}{b2}{b3}'
                if _pattern_valid(pattern, pos_cdns):
                    valid_patterns.append(pattern)

    # reverse back
    if reverse:
        valid_patterns = [vp[::-1] for vp in valid_patterns]

    return valid_patterns


def _patterns_to_overhangs(AA1_patterns, AA2_patterns, acceptable_overhangs):
    """ Converts patterns from the two codon sets into the possible overhangs.

    Essentially finds combinations of patterns from the pattern lists that are
    length 4 when wildcards are eliminated. These combination are valid Golden
    Gate overhangs. The number of left wildcards is recorded with the overhang
    to represent the position of the overhang relative to the two breakpoint
    codons.

    Args:
        AA1_patterns: List of valid patterns for first residue at breakpoint.
            Wildcard '.' should be at the left of patterns.
        AA2_patterns: List of valid patterns for second residue at breakpoint.
            Wildcard '.' should be at the right of patterns.
        acceptable_overhangs: List of overhangs that pass the criteria in
            acceptable_overhang(). Each overhang is of the form
            (overhang index, overhang sequence).

    Returns:
        List of overhangs that are valid for the given breakpoint.
    """
    # bin patterns by length
    AA1_binned = [[], [], []]
    AA2_binned = [[], [], []]
    for pat in AA1_patterns:
        pat = pat.replace('.', '')
        AA1_binned[len(pat) - 1].append(pat)
    for pat in AA2_patterns:
        pat = pat.replace('.', '')
        AA2_binned[len(pat) - 1].append(pat)
    AA1_binned = reversed(AA1_binned)  # reverse for combination with AA2_pats

    # find any length 4 overhangs (e.g. all len 1 AA1 pats w/ len 3 AA pats)
    # pos is 3 - len(AA1_pat)
    overhangs = []
    for pos, (AA1_pats, AA2_pats) in enumerate(zip(AA1_binned, AA2_binned)):
        if not AA1_pats or not AA2_pats:
            continue
        for pat1, pat2 in itertools.product(AA1_pats, AA2_pats):
            oh = pat1 + pat2
            if oh in acceptable_overhangs:
                overhangs.append((pos, oh))
    return overhangs


def find_GG_breakpoints(alignment, vector_overhangs):
    """ Finds valid Golden Gate sites to place breakpoints.

    Valid places to put breakpoints are where Golden Gate assembly can occur
    between adjacent blocks from any two parents. These sites may be considered
    as potential block endpoints for the SCHEMA-RASPP algorithm. Possible
    overhangs for each breakpoint are also returned as dictionary values.

    Args:
        alignment: List of tuples that contain the amino acid of each parent
            in the MSA at the position corresponding to the tuple index.
        vector_overhangs: Tuple of tuples containing information about the
            vector overhangs. First and second elements of outer tuple are the
            start and end overhangs, respectively. Each overhang is of the form
            (overhang index, overhang sequence).

    Returns:
        Dictionary of int: list pairs where keys are valid breakpoint indices
            and values are a list of valid overhangs for the given breakpoint
            represented by (overhang index, overhang sequence) tuples.
    """
    acceptable_overhangs = [oh for oh in all_overhangs
                            if acceptable_overhang(oh, vector_overhangs)]

    breakpoints = {0: [vector_overhangs[0]]}  # start with beginning

    for bp_num in range(1, len(alignment)):
        AA1 = set(alignment[bp_num-1])  # the last postion of block
        AA2 = set(alignment[bp_num])  # the first position of the next block

        if '-' in AA1 or '-' in AA2:  # can't recombine in a loop
            continue

        # codons that could give each AA in AA1 and AA2
        pos1_codons = [general_tools.rev_code[p] for p in AA1]
        pos2_codons = [general_tools.rev_code[p] for p in AA2]

        # overhang patterns that can give all AAs
        pos1_patterns = _get_valid_patterns(pos1_codons, reverse=True)
        pos2_patterns = _get_valid_patterns(pos2_codons)

        bp_overhangs = _patterns_to_overhangs(pos1_patterns, pos2_patterns,
                                              acceptable_overhangs)
        if not bp_overhangs:
            continue

        breakpoints[bp_num] = bp_overhangs

    breakpoints[len(alignment)] = [vector_overhangs[1]]

    return breakpoints


def generate_weighted_E_matrix(alignment, weighted_contacts):
    """ Generate the SCHEMA energy matrix.

    The matrix is a square matrix with the same size as the length of
    alignment. The calculation is weighted by the values in the contacts dict.
    This could be the contact frequency. To get the unweighted behaviour, just
    set all contact weights=1.

    Args:
        alignment: List of tuples that contain the amino acid of each parent
            in the MSA at the position corresponding to the tuple index.
        weighted_contacts: Dictionary of tuple: float pairs where keys are
            a tuple with the indices of the two contacting residues and values
            are the weight associated with the contacts.

    Returns:
        SCHEMA energy matrix, used for calculating the SCHEMA energy of blocks.
    """
    E_matrix = numpy.zeros((len(alignment), len(alignment)))
    for (i, j), weight in weighted_contacts.items():
        AAs_i, AAs_j = alignment[i], alignment[j]
        parental = set(zip(AAs_i, AAs_j))  # get parental pairs of AAs at i,j
        broken = 0
        for combo in itertools.product(AAs_i, AAs_j):
            if combo not in parental:
                broken += 1  # contact is broken, record it
        E_matrix[i, j] = (broken * weight) / (len(AAs_i)**2)
    return E_matrix


def generate_blocks(breakpoints, minBL=10, maxBL=100):
    """Finds all blocks allowed by the block constraints minBL and maxBL.

    Args:
        breakpoints: Dictionary of int: list pairs from find_GG_breakpoints.
        minBL: Int minimum length of blocks.
        maxBl: Int maximum length of blocks.

    Returns:
        List of tuples that contain the start and end breakpoints of each
            block found.
    """
    blocks = set()
    for bp1, bp2 in itertools.combinations(breakpoints, r=2):
        if minBL <= bp2-bp1 <= maxBL:
            blocks.add((bp1, bp2))
    return blocks


def _block_energy(start, end, E_matrix):
    """ Calculate the SCHEMA energy of a block.

    Args:
        start: Int block's start breakpoint index.
        end: Int block's end breakpoint index.
        E_matrix: SCHEMA energy matrix from generate_weighted_E_matrix.
    """
    energy = E_matrix[start:end, end:].sum()
    return energy


# Represents node in networkx graph.
Node = collections.namedtuple('Node', ['col', 'bp'])


def _graph_column(G, blocks, E_matrix, curr_col):
    """ Add a column to the SCHEMA-RASPP graph.

    Column i represents the ith chimera breakpoint. Nodes within the column
    represent the breakpoint position within the sequence. Edges entering
    this column are the ith blocks in the chimeric sequences and are weighted
    by SCHEMA energy. Adding a column represents adding a new block to the
    possible chimeric sequences.

    Args:
        G: SCHEMA-RASPP graph.
        blocks: List of tuples from generate_blocks.
        E_matrix: SCHEMA energy graph from generate_weighted_E_matrix.
        column: Int index of column to be added.
    """
    prev_col = curr_col - 1
    previous_bps = set(node.bp for node in G if node.col == prev_col)
    for start, end in blocks:
        if start in previous_bps:
            start_node = Node(prev_col, start)
            end_node = Node(curr_col, end)
            G.add_edge(start_node, end_node, weight=_block_energy(start, end,
                                                                  E_matrix))


def _remove_deadends(G):
    """ Remove edges that do not have output edges.

    Args:
        G: SCHEMA_RASPP graph.
    """
    num_bl = max(n.col for n in G)
    for col in range(num_bl-1, 0, -1):
        col_bps = set(node.bp for node in G if node.col == col)
        for bp in col_bps:
            node = Node(col, bp)
            if G.out_degree(node) == 0:
                G.remove_node(node)


def _construct_graph(num_bl, blocks, E_matrix):
    """ Make the SCHEMA-RASPP graph.

    Represents the space of chimeric sequence as a directed graph. Graph nodes
    represent breakpoints are arranged in columns represented breakpoint
    indices in the chimeric sequences. Edges connect a given node to a node in
    the next column and represent blocks. Paths through the graph represent
    specific breakpoint positions of a resulting chimeric library. View Fig. 2
    from Endelman JB et al. "Site-directed protein recombination as a
    shortest-path problem." (2004) for a visualization.

    Args:
        num_bl: Int number of blocks in the resulting libraries.
        blocks: List of tuples from generate_blocks.
        E_matrix: SCHEMA energy graph from generate_weighted_E_matrix.

    Returns:
        G: constructed SCHEMA-RASPP graph.
        first_node: The root node at column 0 in G.
        last_node: The final node at column num_bl in G.
    """
    G = networkx.digraph.DiGraph()

    # make column 0
    first_node = Node(0, 0)
    G.add_node(first_node)

    # make middle columns
    for column in range(1, num_bl):
        _graph_column(G, blocks, E_matrix, column)

    # make last column
    column = num_bl
    last_bp = max(b[1] for b in blocks)
    ending_blocks = [b for b in blocks if b[1] == last_bp]
    _graph_column(G, ending_blocks, E_matrix, column)
    last_node = [node for node in G if node.col == column][0]

    # can remove nodes that don't have outgoing edges
    _remove_deadends(G)

    return G, first_node, last_node


def _path_energy(G, path):
    """ Calculate the energy of a path in G. Represents the SCHEMA energy of a
    library. """
    return sum(G.get_edge_data(n1, n2)['weight'] for n1, n2
               in zip(path[:-1], path[1:]))


def _remove_blocks(G, cutoff_len, min_max):
    """ Remove blocks from the graph based on a cutoff length.

    Args:
        G: SCHEMA-RASPP graph.
        cutoff_len: Int block length at which exceeding blocks are removed.
        min_max: 'max' or 'min' stating whether cutoff_len is a maximum or
            minimum cutoff.
    """
    removal_edges = []
    for edge in G.edges:
        node1, node2 = edge
        block_length = node2.bp - node1.bp
        if min_max == 'min' and block_length < cutoff_len:
            removal_edges.append(edge)
        elif min_max == 'max' and block_length > cutoff_len:
            removal_edges.append(edge)
    G.remove_edges_from(removal_edges)

    # can remove nodes that don't have outgoing edges
    _remove_deadends(G)

    # remove nodes with no edges
    removal_nodes = [n for n in G if G.degree(n) == 0]
    G.remove_nodes_from(removal_nodes)


def _find_shortest_path(G, first_node, last_node, libraries):
    """ Find the shortest path through G from first_node to last_node and
    update the SCHEMA energy of the libraries with the culmulative edge energy
    through the path. """
    path = networkx.dijkstra_path(G, first_node, last_node)
    energy = _path_energy(G, path)
    bps = tuple(node.bp for node in path)
    libraries[bps] = {'energy': energy}


def shortest_path_recombination(num_bl, blocks, E_matrix):
    """ Find libraries with minimum energy with block length constraints.

    Scan over all pairs of minimum and maximum block lengths, using a
    SCHEMA-RASPP graph and Dijkstra's algorithm to find the minimum energy
    library for each pair.

    Args:
        num_bl: Int number of blocks in the resulting libraries.
        blocks: List of tuples from generate_blocks.
        E_matrix: SCHEMA energy graph from generate_weighted_E_matrix.

    Returns:
        Dictionary of lib_bps: {'E':' <SCHEMA energy>} pairs where lib_bps is a
            tuple of breakpoint indices for a library.
    """

    masterG, first_node, last_node = _construct_graph(num_bl, blocks, E_matrix)

    align_len = max([end for _, end in blocks])
    minBL = min(end - start for start, end in blocks)
    maxBL = max(end - start for start, end in blocks)
    libraries = {}

    # solve shortest path using dijkstra algorithm
    _find_shortest_path(masterG, first_node, last_node, libraries)

    # scan over all minBL+maxBL combinations
    # scan over min block sizes from minBL to seq_len/num_blocks
    equal_size_blocks = align_len // num_bl
    print(f'Total iterations: {equal_size_blocks-minBL} (later iterations'
          ' are faster)', flush=True)
    for min_cutoff in range(minBL, equal_size_blocks + 1):
        print(f'Current iteration: {min_cutoff-minBL}\r', end='', flush=True)
        G = masterG.copy()
        _remove_blocks(G, min_cutoff, min_max='min')
        # scan over max block sizes from maxBL till it can't solve the SPP
        for max_cutoff in range(maxBL, 0, -1):
            _remove_blocks(G, max_cutoff, min_max='max')
            try:
                _find_shortest_path(G, first_node, last_node, libraries)
            except (networkx.exception.NetworkXNoPath,
                    networkx.exception.NodeNotFound):
                # this happens once there are no paths
                break

    print(flush=True)
    return libraries


def _get_block_alignment(lib_bps, alignment):
    """ Assign block numbers to each element in alignment.

    Args:
        lib_bps: Tuple of ints represent library breakpoint indices. Keys from
            libraries.
        alignment: List of tuples that contain the amino acid of each parent
            in the MSA at the position corresponding to the tuple index.

    Returns:
        List of tuples. Same as alignment but each tuple has the block number
            of the associated residue at position 0.
    """
    block_alignment = []
    start_bps, end_bps = lib_bps[:-1], lib_bps[1:]
    for blnum, (start, end) in enumerate(zip(start_bps, end_bps)):
        extension = [(blnum,) + AAs for AAs in alignment[start:end]]
        block_alignment.extend(extension)
    return block_alignment


def _get_binned_blocks(block_alignment):
    """ Bins parental block sequences by position.

    Args:
        block_alignment: List of tuples returned by _get_block_alignment.

    Returns:
        List of lists of block sequences in bins. Bins are represented by each
            inner list, which are indexed according to the block location of
            block sequences they hold.
    """
    block_bins = {}
    for blk_num, *AAs in block_alignment:
        if blk_num in block_bins:
            block_bins[blk_num].append(AAs)
        else:
            block_bins[blk_num] = [AAs]
    block_bins = [zip(*blk_bin) for blk_bin in block_bins.values()]
    binned_blocks = [[''.join(pb) for pb in blk_bin] for blk_bin in block_bins]
    return binned_blocks


def _chimera_muts(chimera_seq, parent_seq):
    """ Calculates the number of mutations required to convert chimera_seq
    to parent_seq. """
    return sum(1 for cAA, pAA in zip(chimera_seq, parent_seq) if cAA != pAA)


def _calculate_M(block_alignment, sample=False):
    """ Calculate the average mutational diversity in a library.

    Args:
        block_alignment: List of tuples returned by _get_block_alignment.

    Return:
        float average number of mutations from nearest parent in library.
    """
    block_seq, *parent_lists = zip(*block_alignment)
    parents = [''.join(par) for par in parent_lists]
    block_nums = sorted(set(block_seq))
    binned_blocks = _get_binned_blocks(block_alignment)

    totalM = 0
    for ch in itertools.product(*binned_blocks):
        chimera_seq = ''.join(ch)
        totalM += min(_chimera_muts(chimera_seq, ps) for ps in parents)
    M = totalM / len(parents) ** len(block_nums)
    return M


def update_M(libraries, alignment):
    """ Update libraries with average number of mutations from nearest parent.

    Args:
        libraries: Dictionary of lib_bps: {'E':' <SCHEMA energy>} pairs where
            lib_bps is a tuple of breakpoint indices for a library. Created
            by shortest_path_recombination.
        alignment: List of tuples that contain the amino acid of each parent
            in the MSA at the position corresponding to the tuple index.
    """
    progress = general_tools.ProgressOutput(len(libraries))
    for itr, (lib_bp, lib_attrs) in enumerate(libraries.items()):
        progress.update(itr)
        block_alignment = _get_block_alignment(lib_bp, alignment)
        M = _calculate_M(block_alignment)
        lib_attrs['M'] = M


def _calculate_GG_prob(overhang_seqs):
    """ Golden Gate ligation probability.

    Calculates the pseudo-probability of correct GG ligation given a list of
    overhang sequences. Ligation data is from Potapov V et al 2018. This isn't
    a real probability but scales from 0 to 1 and is roughly representative
    of the compatibility of the overhangs in the set.

    Args:
        overhang_seqs: List of str representing overhang DNA sequences.

    Returns:
        float probability of GG ligation.
    """
    complements = [complementary_sequence(oh) for oh in overhang_seqs]

    all_overhang_seqs = overhang_seqs + list(reversed(complements))

    # construct matrix of pairwise overhang ligation counts
    overhang_set_counts = []
    for oh1 in all_overhang_seqs:
        oh_set_counts_row = []
        for oh2 in reversed(all_overhang_seqs):
            oh_set_counts_row.append(normalized_ligation_counts[oh1][oh2])
        overhang_set_counts.append(oh_set_counts_row)

    # find product of matrix diagonals each divided by their row sum
    GG_prob = 1
    for i in range(len(overhang_set_counts)):
        GG_prob *= (overhang_set_counts[i][i] / sum(overhang_set_counts[i]))

    return GG_prob


def update_GG_prob(libraries, breakpoints):
    """ Updates libraries with Golden Gate information.

    Each library is updated with the set of overhangs with highest GG_prob and
    that GG_prob. Libraries with low GG_prob are likely to have issues in
    Golden Gate reactions.

    Args:
        libraries: Dictionary of lib_bps: {'E':' <SCHEMA energy>} pairs where
            lib_bps is a tuple of breakpoint indices for a library. Created
            by shortest_path_recombination.
        breakpoints: Dictionary of int: list pairs from find_GG_breakpoints.
    """
    progress = general_tools.ProgressOutput(len(libraries))
    for itr, (bp_positions, lib_attrs) in enumerate(libraries.items()):
        progress.update(itr)
        lib_breakpoints = {bp_pos: breakpoints[bp_pos] for bp_pos
                           in bp_positions}
        max_gg_prob = 0
        for oh_set in itertools.product(*lib_breakpoints.values()):
            overhangs = [pos_oh[1] for pos_oh in oh_set]
            gg_prob = _calculate_GG_prob(overhangs)
            if gg_prob > max_gg_prob:  # better overhang set found
                max_gg_prob = gg_prob
                max_set = oh_set
                if gg_prob == 1:
                    break

        lib_attrs.update({'GG_prob': max_gg_prob})
        if lib_attrs['GG_prob'] != 0:
            lib_attrs['GG_sites'] = max_set

    return libraries
