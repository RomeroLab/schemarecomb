import pickle
import networkx
import numpy
import itertools
from Bio.Data.CodonTable import standard_dna_table

# Set run parameters.
minBL = 40                                  # Block length restrictions
maxBL = 400
num_bl = 8
alignment_fn = 'alignment.p'
contacts_fn = 'contacts.p'
libraries_out_fn = 'libraries.p'
threshold = 245                             # Minimum threshold for overhangs read from data
vector_overhangs = ['TATG', 'TGAG']
overhang_length = 4

code = standard_dna_table.forward_table.copy()

# Remove rare codons from block junction selection, found from various sources online.
rare = ('ATA', 'AGG', 'TCG', 'ACG', 'AGA', 'GGA', 'CTA', 'GGG', 'CCC', 'CGG', 'CGA') 
(code.pop(r) for r in rare)
rev_code = dict((a,tuple([k for k in code.keys() if code[k]==a])) for a in set(code.values()))
rev_code['-']=()

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def complementary_sequence(sequence):
    return ''.join(reversed([complement[bp] for bp in sequence]))

def acceptable_overhang(overhang, normalized_ligation_counts):
    """Checks that an overhang is valid for GoldenGate cloning, e.g. not complementary, not full-GC,
    above threshold, and is not a vector_overhang.

    Args:
        overhang (str): 5'->3' string of base-pairs that function as sticky ends in GoldenGate cloning.

    Returns:
        bool: True if overhang is valid, false if invalid.
    """
    if overhang == complementary_sequence(overhang):
        return False
    if 'A' not in overhang and 'T' not in overhang:
        return False
    if normalized_ligation_counts[overhang][complementary_sequence(overhang)] < threshold:
        return False
    if overhang in vector_overhangs or complementary_sequence(overhang) in vector_overhangs:
        return False
    return True

def possible_GG_site(breakpoint_num, alignment, normalized_ligation_counts):
    """Collects all valid overhangs at amino acids breakpoint_num-1 and breakpoint_num in alignment by 
    considering all possible codon changes.

    Args:
        breakpoint_num (int): Identifies location of breakpoint to check.  Corresponds to the alignment index 
            of the second amino acid in the breakpoint amino acid pair.
        alignment (list of tuples of strings): Multiple sequence alignment (MSA) formatted as a list of tuples
            where each tuple contains the ordered parental amino acids at the corresponding MSA position.

    Returns:
        list of tuples: Valid overhangs where the 0th tuple element is the start position of the overhang relative
            to the start position of the first breakpoint codon (0, 1, or 3), and the 1st element is the overhang
            base pairs.
"""
    possible_sites = []
    AA1 = set(alignment[breakpoint_num-1]) # the last postion of block
    AA2 = set(alignment[breakpoint_num])   # the first position of the next block
    if '-' not in AA1 and '-' not in AA2:
        # All combinations of codons that could give AAs at pos1 and pos2.
        combo1 = set(itertools.product(*[rev_code[p] for p in AA1]))
        combo2 = set(itertools.product(*[rev_code[p] for p in AA2]))
        for c1 in combo1:
            for c2 in combo2:
                # Find the number of overlapping bases at end of c1 and beginning of c2.
                end = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c1)]).split('0')[-1])
                beg = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c2)]).split('0')[0])
                if end+beg >= overhang_length:
                    string = c1[0]+c2[0]
                    for i in range(3-end,beg):
                        if (i, string[i:i+4]) not in possible_sites and acceptable_overhang(string[i:i+4], normalized_ligation_counts):
                            possible_sites.append((i, string[i:i+4]))
    return possible_sites

def find_GG_breakpoints(alignment, normalized_ligation_counts):
    """Finds all possible breakpoints and the associated possible overhangs.
    
    Args:
        alignment (list of tuples of strings): Multiple sequence alignment (MSA) formatted as a list of tuples
            where each tuple contains the ordered parental amino acids at the corresponding MSA position. 

    Returns:
        dictionary: Keys are valid breakpoints and values are tuples of the associated possible overhangs where 
            the 0th tuple element is the start position of the overhang relative to the start position of the 
            first breakpoint codon (0, 1, or 3), and the 1st element is the overhang base pairs.
    """
    breakpoints = {0: (1, vector_overhangs[0])}
    for breakpoint_num in range(1, len(alignment)):
        print('\rFinding GG breakpoints: {} out of {}'.format(breakpoint_num, len(alignment)), end='', flush=True)
        possible_sites = possible_GG_site(breakpoint_num, alignment, normalized_ligation_counts)
        if len(possible_sites) > 0:
            breakpoints.update({breakpoint_num: possible_sites})
    print()
    breakpoints.update({len(alignment): (0, vector_overhangs[1])})
    return breakpoints

def weighted_E_matrix(alignment, weighted_contacts):
    """Generates the E matrix weighted by the values in the contacts dict.  To get the original (unweighted)
    behaviour, just set all contact weights to 1

    Args:
        alignment (list of tuples of strings): Multiple sequence alignment (MSA) formatted as a list of tuples
            where each tuple contains the ordered parental amino acids at the corresponding MSA position.
        weighted_contacts (dictionary): Keys are length 2 tuples containing indices of contacting pairs of
            residues and values are the associated weight of the contacts.  This could be the contact frequency
            in multiple PDB files.

    Returns:
       2D numpy array: Represents the SCHEMA energy generated by the recombinational breakage of a given contact.
    """
    E_matrix = numpy.zeros((len(alignment), len(alignment)))
    pariter = range(len(alignment[0]))
    for cnt in weighted_contacts:
        parental = set((alignment[cnt[0]][p], alignment[cnt[1]][p]) for p in pariter)
        broken = len([1 for p1 in pariter for p2 in pariter if (alignment[cnt[0]][p1], alignment[cnt[1]][p2]) 
                   not in parental])
        E_matrix[cnt[0], cnt[1]] = (broken * weighted_contacts[cnt]) / float(len(pariter)**2)
    return E_matrix

def generate_blocks(breakpoints, minBL=10, maxBL=100):
    """Finds all blocks that allowed by the block constraints.

    Args:
        breakpoints (dictionary): Keys are valid breakpoints and values are tuples of the associated possible 
            overhangs where the 0th tuple element is the start position of the overhang relative to the start position 
            of the first breakpoint codon (0, 1, or 3), and the 1st element is the overhang base pairs.  (Only keys 
            required for this method.)
        minBL (int): Minimum length of blocks in amino acids.
        maxBL (int): Maximum length of blocks in amino acids. 

    Returns:
        list of tuples: List elements contain the start and end indices of possible blocks. 
    """
    blocks = set()
    for bp1 in breakpoints:
        for bp2 in breakpoints:
            if (bp1 + minBL) <= bp2 and (bp1 + maxBL) >= bp2:
                blocks.add((bp1, bp2))
    return blocks

def shortest_path_recombination(num_bl, blocks, E_matrix, fast=False):
    """Finds libraries: combinations of blocks that split the aligned sequences. The slow version scans over all pairs
    of minBL and maxBL, which is very slow, but more thorough.  The fast version incrementally increases minBL and 
    then incrementally decreases maxBL.

    Args:
        num_bl (int): Number of blocks in a library.
        blocks (list of tuples): List elements contain the start and end indices of possible blocks.
        E_matrix (2D numpy array): Represents the SCHEMA energy generated by the recombinational breakage of a given 
            contact.
        fast (bool): False for slow version, True for fast version.

    Returns:
        dictionary: Keys are tuples of breakpoints.  Values are a dictionary containing the calulated SCHEMA energy.
    """
    print('Entering shortest_path_recombination', flush=True)

    block_len = dict((bl, bl[1]-bl[0]) for bl in blocks)
    block_energy = dict((bl, E_matrix[bl[0]:bl[1], bl[1]:].sum()) for bl in blocks)

    align_len = max([max(b) for b in blocks])
    minBL = min(block_len.values())
    maxBL = max(block_len.values())
    G = networkx.digraph.DiGraph()
    libraries = dict()

    # Build in the first num_bl-1 blocks, node = (bl,(start,end)), and edges are added between compatible nodes: 
    #    (bl,(start,end)) and (bl+1,(end,bl2end))
    G.add_node((0,(0,0)))
    for bl in range(1,num_bl):
        print('\rLoop 1: {} remaining out of {}'.format(bl, num_bl), flush=True, end='')
        previous = set(n for n in G if (n[0]+1)==bl)
        for node in previous:
            G.add_weighted_edges_from([(node, (bl, b), block_energy[b]) for b in blocks if node[1][1]==b[0]])

    # Add in last block so that it ends exactly at align_len.
    print('\nLoop 2', flush=True)
    bl += 1 
    previous = set(n for n in G if (n[0]+1)==bl)
    for node in previous:
        if (node[1][1], align_len) in blocks:
            G.add_edge(node, (bl, (node[1][1], align_len)), weight=block_energy[(node[1][1], align_len)])
            G.add_edge((bl, (node[1][1], align_len)),(bl+1, (align_len, align_len)), weight=0)

    # Solve shortest path using dijkstra algorithm.
    path = networkx.dijkstra_path(G, (0, (0,0)), (num_bl+1, (align_len, align_len)))
    energy = sum([G.get_edge_data(path[i], path[i+1])['weight'] for i in range(len(path)-1)])
    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy': energy}

    masterG = G.copy()

    if fast:
        # Make the shortest allowed block incrementally larger.
        for sub in range(10000):
            print('\rLoop 3: {} remaining out of {}'.format(sub, 10000), flush=True, end='')
            remove = set(b for b in blocks if block_len[b]==minBL+sub)
            G.remove_nodes_from([n for n in G if n[1] in remove])
            try:
                path = networkx.dijkstra_path(G, (0, (0,0)), (num_bl+1, (align_len, align_len)))
                energy = sum([G.get_edge_data(path[i], path[i+1])['weight'] for i in range(len(path)-1)])
                libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
            except: # this happens once there are no paths
                break
        print()

        G = masterG
        # make the largest allowed block incrementally smaller
        for sub in range(10000):
            print('\rLoop 4: {} remaining out of {}'.format(sub, 10000), flush=True, end='')
            remove = set(b for b in blocks if block_len[b]==maxBL-sub)
            G.remove_nodes_from([n for n in G if n[1] in remove])
            try:
                path = networkx.dijkstra_path(G, (0, (0,0)), (num_bl+1, (align_len, align_len)))
                energy = sum([G.get_edge_data(path[i], path[i+1])['weight'] for i in range(len(path)-1)])
                libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy': energy}
            except: # this happens once there are no paths
                break
        print()

    else:
	# Scan over all minBL+maxBL combinations.
        for min_sub in range(int((align_len/num_bl))-minBL):
            print('\rLoop 4: {}} remaining out of {}'.format(min_sub, int((align_len/num_bl))-minBL), 
                     flush=True, end='')
            G = masterG.copy()
            remove = set(b for b in blocks if block_len[b]<=(minBL+min_sub))
            G.remove_nodes_from([n for n in G if n[1] in remove])
            for max_sub in range(10000): # Scan over max block sizes from maxBL till it can't solve the SPP.
                remove = set(b for b in blocks if block_len[b]>=(maxBL-max_sub))
                G.remove_nodes_from([n for n in G if n[1] in remove])
                try:
                    path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
                    energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
                    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
                except:       # this happens once there are no paths
                    break
    return libraries

def calculate_M(alignment, sample=False):
    """Calculates the expected number of mutations in each chimera relative to the closest parent.
     
    Args:
        alignment (list of tuples of strings): Multiple sequence alignment (MSA) formatted as a list of tuples
            where each tuple contains the ordered parental amino acids at the corresponding MSA position.
        sample (False or int): If False, the number of mutations for every possible chimera is counted in order to 
            calculate M.  If int, numbers of mutations for that many randomly chosen chimeras are counted.

    Returns:
        int: (Calculated) Expected number of mutations in each chimera relative to the closest parent.
    """
    block_alignment = []
    for i in range(len(library)-1):
        block_alignment.extend([(i+1,)+tuple(p) for p in alignment[library[i]:library[i+1]]])

    parents = [''.join(s) for s in tuple(zip(*block_alignment))[1:]]
    blocks = sorted(set([p[0] for p in block_alignment]))
    totalM = 0
    if not sample:
        all_chimeras = itertools.product(*[[''.join(s) for s in zip(*[b[1:] for b in block_alignment if b[0]==bl])] 
                for bl in blocks])
        for ch in all_chimeras:
            totalM += min([sum([p[i]!=''.join(ch)[i] for i in range(len(block_alignment))]) for p in parents])
        M = float(totalM)/((len(block_alignment[0])-1)**len(blocks))
    else:
        for rnd in range(sample):
            ch = [choice([''.join(s) for s in zip(*[b[1:] for b in block_alignment if b[0]==bl])]) for bl in blocks]
            totalM += min([sum([p[i]!=''.join(ch)[i] for i in range(len(block_alignment))]) for p in parents])
        M = float(totalM)/sample
    return M

def update_M(libraries, alignment, sample=False):
    """Adds the M values to the libraries dict.

    Args:
        libraries (dict): Keys are a tuple of breakpoints for the library. Values are dictionaries that include the
            SCHEMA energy value calculated for the library.
        alignment (list of tuples of strings): Multiple sequence alignment (MSA) formatted as a list of tuples
            where each tuple contains the ordered parental amino acids at the corresponding MSA position.
        sample (False or int): If False, the number of mutations for every possible chimera is counted in order to 
            calculate M.  If int, numbers of mutations for that many randomly chosen chimeras are counted.

    Returns:
        dict: libraries arg with values dictionaries updated with the calculated M value.
    """
    print('Updating M', flush=True)
    for lib in libraries:
        M = calculate_M(library_to_block_alignment(lib,alignment),sample)
        libraries[lib]['M'] = M
    return libraries

def calculate_GG_prob(list_of_overhangs, normalized_ligation_counts):
    """Estimates the likelihood of correct GG ligation given a set of overhangs using the overhang ligation data from 
    Lohman et al. 2018.  This estimate is not a direct calculation of ligation probability.


    Args:
        list_of_overhangs (list): Set of overhangs on one DNA strand. Complementary overhangs are added automatically.
        normalized_ligation_counts (dict): 2D dictionary of all Lohman et al. 2018 ligation counts normalized 100,000
            total observations.  Keys in both levels of dict are overhang strings while final values are floats
            representing normalized ligation counts between the overhangs in the keys.

    Returns:
        float: The "overhang probability" is calculated for each overhang by dividing the number of observed correct 
            ligations by the total number of ligations involving that overhang.  The "probabilities" for each overhang
            are then multiplied to get the overhang set GG_prob.
    """
    complements = []
    for oh in list_of_overhangs:
        complements.append(complementary_sequence(oh))
    all_overhangs = list_of_overhangs + list(reversed(complements))
    overhang_set_counts = []
    for oh1 in all_overhangs:
        overhang_set_counts_row = []
        for oh2 in reversed(all_overhangs):
            overhang_set_counts_row.append(normalized_ligation_counts[oh1][oh2])
        overhang_set_counts.append(overhang_set_counts_row)
    GG_prob = 1
    for i in range(len(overhang_set_counts)):
        GG_prob *= (overhang_set_counts[i][i] / sum(overhang_set_counts[i]))
    return GG_prob
    
def update_GG_prob(libraries, alignment, breakpoints, normalized_ligation_counts):
    """Updates libraries with the set of overhangs with highest GG_prob.
    
    Args:
        libraries (dict): Keys are a tuple of breakpoints for the library. Values are dictionaries that include the
            SCHEMA energy value and M calculated for the library.
        alignment (list of tuples of strings): Multiple sequence alignment (MSA) formatted as a list of tuples
            where each tuple contains the ordered parental amino acids at the corresponding MSA position.
        breakpoints (dictionary): Keys are valid breakpoints and values are tuples of the associated possible 
            overhangs where the 0th tuple element is the start position of the overhang relative to the start position 
            of the first breakpoint codon (0, 1, or 3), and the 1st element is the overhang base pairs.
        normalized_ligation_counts (dict): 2D dictionary of all Lohman et al. 2018 ligation counts normalized 100,000
            total observations.  Keys in both levels of dict are overhang strings while final values are floats
            representing normalized ligation counts between the overhangs in the keys.

    Returns:
        dict: libraries arg with values dictionaries updated with the calculated GG_prob value and best set of 
            overhangs.
    """
    all_possible_sites = {}
    for itr, lib in enumerate(libraries):              
        lib_valid = True
        lib_possible_sites = {}
        for breakpoint_num in lib[1:-1]:
            if len(all_possible_sites[breakpoint_num]) == 0:
                lib_valid = False
            lib_possible_sites.update({breakpoint_num: breakpoints[breakpoint_num]})
        if not lib_valid:
            libraries[lib].update({'GG_prob': 0})
            print('not valid')
            continue
        max_p = 0
        max_set = ()
        prod = list(itertools.product(*lib_possible_sites.values()))
        print('%i out of %i libraries: %i possible GG sites combinations' % (itr, len(libraries), len(prod)), 
                flush=True)
        for s in prod:
            p = calculate_GG_prob(list(list(zip(*s))[1]) + vector_overhangs, normalized_ligation_counts)
            if p > max_p:
                max_p = p
                max_set = s
                if p == 1:
                    break
        libraries[lib].update({'GG_prob':max_p})
        if libraries[lib]['GG_prob'] != 0:
            libraries[lib].update({'GG_sites':max_set})
    print("Done updating GG_prob")        
    return libraries

if __name__ == '__main__':
        normalized_ligation_counts = pickle.load(open('normalized_ligation_counts_18h_37C.p', 'rb'))
        alignment = pickle.load(open(alignment_fn, 'rb'))
        weighted_contacts = pickle.load(open(contacts_fn, 'rb'))

	# A breakpoint specifies the first position of a new block.
        breakpoints = find_GG_breakpoints(alignment, normalized_ligation_counts)

	# The same E, but weighted by the values in the contacts dict. This could be the contact frequency. 
        # To get unweighted, just set all contact weight==1""
        E_matrix = weighted_E_matrix(alignment, weighted_contacts)

	# generate all allowed blocks
        blocks = generate_blocks(breakpoints, minBL, maxBL)

	# run RASPP
        libraries = shortest_path_recombination(num_bl, blocks, E_matrix, fast=False)

	# Add M values to the libraries dictionary
        update_M(libraries, alignment)

        pickle.dump(libraries, open(libraries_fn, 'wb'))

	#Requires 'normalized_ligation_counts_18h_37C.p' file to be in same dir
        libraries = update_GG_prob(libraries, alignment)

        pickle.dump(libraries, open('gg_' + libraries_fn, 'wb'))
