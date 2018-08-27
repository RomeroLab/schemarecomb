import pickle
import networkx
import numpy
import itertools

# define the library properties min block length, max block length, number of blocks, and parents
minBL = 40
maxBL = 400
num_bl = 8

alignment_fn = 'alignment.p'
contacts_fn = 'contacts.p'
libraries_fn = 'libraries.p'

alignment = pickle.load(open(alignment_fn, 'rb'))
contacts = pickle.load(open(contacts_fn, 'rb'))
normalized_ligation_counts = pickle.load(open(ligation_counts_fn,'rb'))

threshold = 245
vector_overhangs = ['TATG','TGAG']
overlap=4

code = {'AAA': 'K','AAC': 'N','AAG': 'K','AAT': 'N','ACA': 'T','ACC': 'T','ACG': 'T','ACT': 'T','AGA': 'R','AGC': 'S','AGG': 'R','AGT': 'S','ATA': 'I','ATC': 'I','ATG': 'M','ATT': 'I',
        'CAA': 'Q','CAC': 'H','CAG': 'Q','CAT': 'H','CCA': 'P','CCC': 'P','CCG': 'P','CCT': 'P','CGA': 'R','CGC': 'R','CGG': 'R','CGT': 'R','CTA': 'L','CTC': 'L','CTG': 'L','CTT': 'L',
        'GAA': 'E','GAC': 'D','GAG': 'E','GAT': 'D','GCA': 'A','GCC': 'A','GCG': 'A','GCT': 'A','GGA': 'G','GGC': 'G','GGG': 'G','GGT': 'G','GTA': 'V','GTC': 'V','GTG': 'V','GTT': 'V',
        'TAA': '*','TAC': 'Y','TAG': '*','TAT': 'Y','TCA': 'S','TCC': 'S','TCG': 'S','TCT': 'S','TGA': '*','TGC': 'C','TGG': 'W','TGT': 'C','TTA': 'L','TTC': 'F','TTG': 'L','TTT': 'F',}

# remove rare codons--we don't want to allow these at the block junctions (I found these from various sources online)
rare = ('ATA', 'AGG', 'TCG', 'ACG', 'AGA', 'GGA', 'CTA', 'GGG', 'CCC', 'CGG', 'CGA') 
[code.pop(r) for r in rare]
aminos = set(code.values())
rev_code = dict((a,tuple([k for k in code.keys() if code[k]==a])) for a in aminos)
rev_code['-']=()

complement = {'A':'T','C':'G','G':'C','T':'A'}

def complementary_sequence(sequence):
    return ''.join(reversed([complement[bp] for bp in sequence]))

def acceptable_overhang(overhang):
    """Returns true if overhang is not complementary, full-GC, below threshold, or is not a vector_overhang."""
    if overhang == complementary_sequence(overhang):
        return False
    if 'A' not in overhang and 'T' not in overhang:
        return False
    if normalized_ligation_counts[overhang][complementary_sequence(overhang)] < threshold:
        return False
    if overhang in vector_overhangs or complementary_sequence(overhang) in vector_overhangs:
        return False
    return True

def breakpoint_good(breakpoint_num):
    AA1 = set(alignment[breakpoint_num-1]) # the last postion of block
    AA2 = set(alignment[breakpoint_num]) # the first position of the next block
    combo1 = set(itertools.product(*[rev_code[p] for p in AA1])) # all combinations of codons that could give AAs at pos1
    combo2 = set(itertools.product(*[rev_code[p] for p in AA2])) # all combinations of codons that could give AAs at pos2
    for c1 in combo1:
        for c2 in combo2:
            end = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c1)]).split('0')[-1]) # number of overlapping bases at end of c1
            beg = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c2)]).split('0')[0]) # number of overlapping bases at beginning of c2
            if end+beg >= overlap:
                string = c1[0]+c2[0]
                for i in range(3-end,beg):
                    if acceptable_overhang(string[i:i+4]):
                        return True
    return False

def find_GG_breakpoints(alignment):
    """Marks breakpoints as possible GG sites if it has a possible GG overhang."""
    breakpoints = [0] # start with beg and end
    for breakpoint_num in range(1, len(alignment)):
        breakpoint_good = False
        AA1 = set(alignment[breakpoint_num-1]) # the last postion of block
        AA2 = set(alignment[breakpoint_num]) # the first position of the next block
        if '-' not in AA1 and '-' not in AA2 and breakpoint_good(breakpoint_num):
            breakpoints.append(breakpoint_num)
    breakpoints.append(len(alignment))
    return breakpoints

def generate_weighted_E_matrix(alignment, weighted_contacts):
    """this is the same as E, but weighted by the values in the contacts dict.  This could be the contact frequency
    To get the original (unweighted) behaviour, just set all contact weights==1"""
    E_matrix = numpy.zeros((len(alignment), len(alignment)))
    pariter = range(len(alignment[0]))
    for cnt in weighted_contacts:
        parental = set((alignment[cnt[0]][p], alignment[cnt[1]][p]) for p in pariter)
        broken = len([1 for p1 in pariter for p2 in pariter if (alignment[cnt[0]][p1], alignment[cnt[1]][p2]) not in parental])
        E_matrix[cnt[0],cnt[1]] = (broken * weighted_contacts[cnt]) / float(len(pariter)**2)
    return E_matrix

def generate_blocks(break_points,minBL=10,maxBL=100):
    """finds all blocks that allowed by the block constraints minBL and maxBL"""
    blocks = set()
    for bp1 in break_points:
        for bp2 in break_points:
            if (bp1+minBL)<=bp2 and (bp1+maxBL)>=bp2:
                blocks.add((bp1,bp2))
    return blocks

def shortest_path_recombination(num_bl, blocks, E_matrix, fast=True):
    """this version scans over all pairs of minBL and maxBL.  This is very slow, but more thorough.  At least the memory requirements are small: if you give it time it will finish"""
    print('Entering shortest_path_recombination',flush=True)

    block_len = dict((bl,bl[1]-bl[0]) for bl in blocks)
    block_energy = dict((bl,E_matrix[bl[0]:bl[1],bl[1]:].sum()) for bl in blocks)

    align_len = max([max(b) for b in blocks])
    minBL = min(block_len.values())
    maxBL = max(block_len.values())
    if verbose: print('filling recombination graph ....')
    G = networkx.digraph.DiGraph()
    libraries = dict()

    print('Loop 1', flush=True)

    # build in the first num_bl-1 blocks, node = (bl,(start,end)), and edges are added between compatible nodes: (bl,(start,end)) and (bl+1,(end,bl2end))
    G.add_node((0,(0,0)))
    for bl in range(1,num_bl):
#        print('\r%i remaining out of %i      ' % (bl, num_bl), flush=True, end='')
        previous = set(n for n in G if (n[0]+1)==bl)
        for node in previous:
            G.add_weighted_edges_from([(node,(bl,b),block_energy[b]) for b in blocks if node[1][1]==b[0]])

    print('\nLoop 2', flush=True)

    # add in last block so that it ends exactly at align_len
    bl += 1 
    previous = set(n for n in G if (n[0]+1)==bl)
    for node in previous:
        if (node[1][1],align_len) in blocks:
            G.add_edge(node,(bl,(node[1][1],align_len)), weight=block_energy[(node[1][1],align_len)])
            G.add_edge((bl,(node[1][1],align_len)),(bl+1,(align_len,align_len)),weight=0)

    # solve shortest path using dijkstra algorithm
    path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
    energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
    if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)

    #fast vs normal difference starts here

    print('Loop 3', flush=True)
     
    # scan over all minBL+maxBL combinations
    masterG = G.copy()

    if fast:
             # make the shortest allowed block incrementally larger
        for sub in range(10000):
            print('\r%i remaining out of %i      ' % (sub, 10000), flush=True, end='')
            remove = set(b for b in blocks if block_len[b]==minBL+sub)
            G.remove_nodes_from([n for n in G if n[1] in remove])
            try:
                path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
                energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
                libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
            except: # this happens once there are no paths
                break

            print('\nLoop 4', flush=True)

            G = masterG
            # make the largest allowed block incrementally smaller
            for sub in range(10000):
                print('\r%i remaining out of %i      ' % (sub, 10000), flush=True, end='')
                remove = set(b for b in blocks if block_len[b]==maxBL-sub)
                G.remove_nodes_from([n for n in G if n[1] in remove])
                try:
                    path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
                    energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
                    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
                except: # this happens once there are no paths
                    break
            print()

    else:
	# scan over all minBL+maxBL combinations
        for min_sub in range(int((align_len/num_bl))-minBL): # scan over min block sizes from minBL to seq_len/num_blocks
	#        print('\r%i remaining out of %i      ' % (min_sub, int((align_len/num_bl))-minBL), flush=True, end='')
            G = masterG.copy()
            remove = set(b for b in blocks if block_len[b]<=(minBL+min_sub))
            G.remove_nodes_from([n for n in G if n[1] in remove])
            for max_sub in range(10000): # scan over max block sizes from maxBL till it can't solve the SPP
	#            if min_sub > int((align_len/num_bl))-minBL - 2:
	#                print('\rmax_sub = %i' % max_sub,flush=True,end='')
                remove = set(b for b in blocks if block_len[b]>=(maxBL-max_sub))
                G.remove_nodes_from([n for n in G if n[1] in remove])
                try:
                    path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
                    energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
                    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
                    if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)
                except: # this happens once there are no paths
                    break
    return libraries

def calculate_M(block_alignment, sample=False):
    """calculates M the normal way, if sample is set to a number, then M is calulated with 'sample' randomly generated chimeras"""
    parents = [''.join(s) for s in tuple(zip(*block_alignment))[1:]]
    blocks = sorted(set([p[0] for p in block_alignment]))
    totalM = 0
    if not sample:
        all_chimeras = itertools.product(*[[''.join(s) for s in zip(*[b[1:] for b in block_alignment if b[0]==bl])] for bl in blocks])
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
    """Adds the M values to the libraries dict.  If sample is a number, then n random samples are used to calculate M """
    print('Updating M', flush=True)
    for lib in libraries:
        M = calculate_M(library2blockalign(lib,alignment),sample)
        libraries[lib]['M'] = M
    return libraries

def possible_GG_site(breakpoint_num, alignment): 
    """Determines whether there is a possible GG overlap at breakpoint_num."""
    possible_sites = []
    AA1 = set(alignment[breakpoint_num-1]) # the last postion of block
    AA2 = set(alignment[breakpoint_num]) # the first position of the next block
    if '-' not in AA1 and '-' not in AA2: # can't recombine in a loop
        combo1 = set(itertools.product(*[rev_code[p] for p in AA1])) # all combinations of codons that could give AAs at pos1
        combo2 = set(itertools.product(*[rev_code[p] for p in AA2])) # all combinations of codons that could give AAs at pos2
        for c1 in combo1:
            for c2 in combo2:
                end = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c1)]).split('0')[-1]) # number of overlapping bases at end of c1
                beg = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c2)]).split('0')[0]) # number of overlapping bases at beginning of c2
                if end+beg >= overlap:
                    string = c1[0]+c2[0]
                    for i in range(3-end,beg):
                        if (i, string[i:i+4]) not in possible_sites and acceptable_overhang(string[i:i+4]):
                            possible_sites.append((i,string[i:i+4]))
    return possible_sites

def calculate_GG_prob(list_of_overhangs):
    """Calculates the probability of correct GG ligation given a list of overhangs."""
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
    
def update_GG_prob(libraries, alignment, vector_overhangs=['TATG','TGAG'], threshold=245):
    """Updates libraries with the set of overhangs with highest GG_prob."""
    all_possible_sites = {}
    for breakpoint_num in range(1,len(alignment)):
        all_possible_sites.update({breakpoint_num : possible_GG_site(breakpoint_num, alignment, vector_overhangs, threshold)})
    for itr, lib in enumerate(libraries):              
        lib_valid = True
        lib_possible_sites = {}
        for breakpoint_num in lib[1:-1]:
            if len(all_possible_sites[breakpoint_num]) == 0:
                lib_valid = False
            lib_possible_sites.update({breakpoint_num : all_possible_sites[breakpoint_num]})
        if not lib_valid:
            libraries[lib].update({'GG_prob':0})
            print('not valid')
            continue
        max_p = 0
        max_set = ()
        prod = list(itertools.product(*lib_possible_sites.values()))
        print('%i out of %i libraries: %i possible GG sites combinations' % (itr, len(libraries), len(prod)), flush=True)
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
	# A breakpoint specifies the first position of a new block.
	breakpoints = find_GG_breakpoints(alignment)

	# The same E, but weighted by the values in the contacts dict. This could be the contact frequency. To get unweighted, just set all contact weight==1""
	E_matrix = generate_weighted_E_matrix(alignment, contacts)

	# generate all allowed blocks
	blocks = generate_blocks(breakpoints, minBL, maxBL)

	# run RASPP
	libraries = shortest_path_recombination(num_bl, blocks, E_matrix, fast=False)

	# Add M values to the libraries dictionary
	update_M(libraries,alignment)

	pickle.dump(libraries, open(libraries_fn, 'wb'))

	#Requires 'normalized_ligation_counts_18h_37C.p' file to be in same dir
	libraries = update_GG_prob(libraries, alignment)

	pickle.dump(libraries, open('gg_' + libraries_fn, 'wb'))
