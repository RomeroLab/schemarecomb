""" Golden Gate SCHEMA-RASPP algorithm.

Run syntax: "python step2.py <num_bl> <minBL> <maxBL> <alignment_fn>
    <contacts_fn> <libraries_fn>"
Example: "python step2.py 8 50 300 alignment.fasta contact.json libraries.json"

Command Line Args:
    num_bl: number of blocks in libraries
    minBL: minimum block length in libraries
    maxBL: maximum blocki length in libraries
    alignment_fn: name of alignment file from Step1
    contacts_fn: name of contacts file from Step1
    libraries_fn: output filename

Output:
    libraries dictionary in library.json

You can also replace all instances of "sys.argv" with the input directly, then
run "python step2.py".

It's recommended to run this script with "nohup <command> &> log.txt &" because
it may take awhile. You can then leave the server and the script will keep
running. Output will be redirected to log.txt.
"""

import json
import sys

from Bio import SeqIO

import raspp_tools

# Read inputs from sys.argv. Can also replace sys.argv below with strings/ints.
num_bl = int(sys.argv[1])  # number of blocks in libraries
minBL = int(sys.argv[2])  # minimum block length in libraries
maxBL = int(sys.argv[3])  # maximum blocki length in libraries
alignment_fn = sys.argv[4]  # name of alignment file from Step1
contacts_fn = sys.argv[5]  # name of contacts file from Step1
libraries_fn = sys.argv[6]  # output filename

print(f'Running SCHEMA-RASPP with the following parameters:', flush=True)
print(f'{num_bl} blocks, minimum block length: {minBL}, maximum block '
      f'length: {maxBL}', flush=True)
print(f'Loaded alignment file: "{alignment_fn}"')
print(f'Loaded contacts file: "{contacts_fn}"')

# Load muscle output fasta file. AA_alignment is a list of tuples where each
# tuple contains the residues at the tuple index position in the MSA.
alignment_SRs = list(SeqIO.parse(alignment_fn, 'fasta'))
AA_seqs = [str(i.seq) for i in alignment_SRs]
AA_alignment = list(zip(*AA_seqs))  # list of AAs at each position

# Load contact dictionary. Dictionary values are 1 to correspond to unweighted
# contact behavior. Change dictionary values to weight contacts.
with open(contacts_fn, 'r') as f:
    contacts = json.load(f)
contacts = {tuple(c): 1 for c in contacts}

# A breakpoint specifies the first position of a new block. breakpoints is a
# dictionary with <breakpoint>: <potential GG overhangs> pairs.
breakpoints = raspp_tools.find_GG_breakpoints(AA_alignment)
print(f'{len(breakpoints)} breakpoints found')

# Generate matrix used to calulcate SCHEMA energy.
E_matrix = raspp_tools.generate_weighted_E_matrix(AA_alignment, contacts)
print('E matrix generated')

# Generate all allowed blocks. blocks is a list of tuples with starting and
# ending breakpoints.
blocks = raspp_tools.generate_blocks(breakpoints, minBL, maxBL)
print('Blocks generated')

# Generate libraries with the minimum value for a given minimum and maximum
# block length combination. libraries is a dictionary of
# lib_bps: {'E':' <SCHEMA energy>} pairs where lib_bps is a tuple of
# breakpoint indices.
print('Starting shortest path recombination to calculate E', flush=True)
libraries = raspp_tools.shortest_path_recombination(num_bl, blocks, E_matrix)
lib_len = len(libraries)
print(f'Shortest path recombination complete: {lib_len} libraries found')

# Add 'M': <average chimera mutations> pairs to the libraries values.
print('Updating M', flush=True)
raspp_tools.update_M(libraries, AA_alignment)
print('M updated', flush=True)

# Add 'GG_prob': <Golden Gate probability> and 'GG_sites': <valid overhangs
# for each lib_bp> pairs the libraries values.
# Requires 'normalized_ligation_counts_18h_37C.p' file to be in code directory.
print('Updating GG probabilities')
libraries = raspp_tools.update_GG_prob(libraries, breakpoints)
print('GG probabilities updated')

# Convert libraries to be json-compatible and save.
libraries = [(k, v) for k, v in libraries.items()]
with open(libraries_fn, 'w') as f:
    json.dump(libraries, f)

print(f'Saved libraries to "{libraries_fn}"')
