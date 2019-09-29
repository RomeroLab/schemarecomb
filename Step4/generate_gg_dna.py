""" Generate DNA fragments from chosen library for ordering.

Run syntax: "python step4.py <AA_alignment_fn> <CDN_alignment_fn>
    <chosen_lib_fn> <frags_order_fn>"
Example: "python step4.py bgl3_AA.fasta bgl3_CDN.fasta bgl3_chosen_lib.json"

Command Line Args:
    AA_alignment_fn: name of aligned amino acid sequences file, from step 1
    CDN_alignment_fn: name of aligned codon sequences file, from step1
    chosen_lib_fn: name of chosen library file, from step 3
    frags_order_fn: name of frags order fasta file

Output:
    DNA fragments for ordering in <frags_order_fn>

You can also replace all instances of "sys.argv" in the code with the input
filenames directly, then run "python step1.py".
"""

import json
import sequence_tools

alignment = pickle.load(open('alignment.p', 'rb'))
breakpoints, lib = pickle.load(open('chosen_lib.p', 'rb'))
vector_overhangs = ['TATG','TGAG']

overhangs = [(2, vector_overhangs[0])] + list(lib['GG_sites']) + [(2, vector_overhangs[1])]

NT_seqs = sequence_tools.read_fasta('DXS_parent_CDN_alignment.fasta')[1]

frags = []

bsa1_start = 'GGTCTC' + 'C'
bsa1_end = 'C' + 'GAGACC'

#library 1
for bp in range(1, len(overhangs) - 1):
    
    start_ind = breakpoints[bp-1]*3 + overhangs[bp-1][0] - 3
    end_ind = breakpoints[bp]*3 + overhangs[bp][0] + 1
    
    new_frags = []
    for seq in NT_seqs:
        raw_frag = seq[start_ind + 4 : end_ind - 4]
        frag_w_overhangs = overhangs[bp-1][1] + raw_frag + overhangs[bp][1]
        frag_w_gg_sites = bsa1_start + frag_w_overhangs + bsa1_end
        new_frags.append(frag_w_gg_sites.replace('-',''))
        
    frags.append(new_frags)
    
    new_frag = overhangs[bp-1][1] + NT_seqs[0][start_ind + 4 : end_ind - 4] + overhangs[bp][1]

print(frags)
