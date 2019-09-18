import json

from Bio import SeqIO

import raspp_tools

# define the library properties min block length, max block length, number of
# blocks, and parents
minBL = 100
maxBL = 400
num_bl = 4

alignment_fn = 'alignment.p'    # name of Step1 output
contacts_fn = 'contacts.p'      # name of Step2 output
libraries_fn = 'libraries.p'    # output filename

# Load muscle output fasta file
alignment_SRs = list(SeqIO.parse(alignment_fn, 'fasta'))
names = [str(i.id) for i in alignment_SRs]
AA_seqs = [str(i.seq) for i in alignment_SRs]
AA_alignment = list(zip(*AA_seqs))  # list of AAs at each position

with open(contacts_fn, 'r') as f:
    contacts = json.load(f)

# A breakpoint specifies the first position of a new block.
breakpoints = raspp_tools.find_GG_breakpoints(AA_alignment)

# The same E, but weighted by the values in the contacts dict. This could be
# the contact frequency. To get unweighted, just set all contact weight==1""
E_matrix = raspp_tools.generate_weighted_E_matrix(AA_alignment, contacts)

# generate all allowed blocks
blocks = raspp_tools.generate_blocks(breakpoints, minBL, maxBL)

# run RASPP
# fast version:
# libraries = raspp_tools.fast_shortest_path_recombination(num_bl, blocks,
#                                                          E_matrix, False)
# slow, but thorough:
libraries = raspp_tools.shortest_path_recombination(num_bl,
                                                    blocks, E_matrix, False)

print('\n Updating M', flush=True)

# Add M values to the libraries dictionary
raspp_tools.update_M(libraries, AA_alignment)

with open(libraries_fn, 'w') as f:
    json.dump(libraries, f)

# Requires 'normalized_ligation_counts_18h_37C.p' file to be in same dir
libraries = raspp_tools.update_GG_prob(libraries, AA_alignment)

with open('gg_' + libraries_fn, 'w') as f:
    json.dump(libraries, f)
