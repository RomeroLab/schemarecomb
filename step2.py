""" Golden Gate SCHEMA-RASPP algorithm.

Run syntax: "python step2.py <num_bl> <minBL> <maxBL> <alignment_fn>
    <contacts_fn> <libraries_fn>"
Example: "python step2.py 8 50 300 bgl3_AA.fasta bgl3_contacts.json
    bgl3_libraries.json"

Command Line Args:
    num_bl: number of blocks in libraries
    minBL: minimum block length in libraries
    maxBL: maximum blocki length in libraries
    start_oh_seq: sequence of start vector overhang
    temp_start_oh_pos: position of start vector overhang relative to first base
        of parent sequence. Can range from -4 to -1, inclusive. Must be
        converted to position relative to imaginary codon at position -1 (i.e.
        add 3) to stay consistent with breakpoint/overhang positioning
        convention
    end_oh_seq: sequence of end vector overhang
    end_oh_pos: position of end vector overhang first base pair relative to
        last codon in parent sequences. Can range from 0 to 3, inclusive
    alignment_fn: name of amino acid alignment file, from step 1
    contacts_fn: name of contacts file, from step 1
    libraries_fn: output filename

Output:
    libraries dictionary in <libraries_fn>, used in step 3

You can also replace all instances of "sys.argv" in the code with the input
filenames directly, then run "python step2.py".

This script generates the graph from Endelman JB et al. "Site-directed protein
recombination as a shortest-path problem." (2004), but uses Dijkstra's
algorithm to solve the shortest-path problem. While this is technically less
efficient, the reduction in number of breakpoints caused by the Golden Gate
criteria makes enumeration of nearly all paths tractable.

It's recommended to run this script with "nohup <command> &> log.txt &" because
it may take awhile. You can then leave the server and the script will keep
running. Output will be redirected to log.txt.

BioPython is a required package. The Romero Lab group server has it installed.
"""

import json
import sys

from Bio import SeqIO

from tools import step2_tools

# Read inputs from sys.argv. Can also replace sys.argv below with strings/ints.
num_bl = int(sys.argv[1])  # number of blocks in libraries
minBL = int(sys.argv[2])  # minimum block length in libraries
maxBL = int(sys.argv[3])  # maximum blocki length in libraries
start_oh_seq = sys.argv[4]  # sequence of start vector overhang
temp_start_oh_pos = int(sys.argv[5])  # position of start vector overhang,
# (must be converted to position relative to 1 imaginary codon back)
start_oh_pos = temp_start_oh_pos + 3
end_oh_seq = sys.argv[6]  # sequence of end vector overhang
end_oh_pos = int(sys.argv[7])  # position of end vector overhang
alignment_fn = sys.argv[8]  # name of alignment file from Step1
contacts_fn = sys.argv[9]  # name of contacts file from Step1
libraries_fn = sys.argv[10]  # output filename

vector_overhangs = ((start_oh_pos, start_oh_seq),
                    (end_oh_pos, end_oh_seq))

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
breakpoints = step2_tools.find_GG_breakpoints(AA_alignment, vector_overhangs)
print(f'{len(breakpoints)} breakpoints found')

# Generate matrix used to calulcate SCHEMA energy.
E_matrix = step2_tools.generate_weighted_E_matrix(AA_alignment, contacts)
print('E matrix generated')

# Generate all allowed blocks. blocks is a list of tuples with starting and
# ending breakpoints.
blocks = step2_tools.generate_blocks(breakpoints, minBL, maxBL)
print('Blocks generated')

# Generate libraries with the minimum value for a given minimum and maximum
# block length combination. libraries is a dictionary of
# lib_bps: {'E':' <SCHEMA energy>} pairs where lib_bps is a tuple of
# breakpoint indices.
print('Starting shortest path recombination to calculate E', flush=True)
libraries = step2_tools.shortest_path_recombination(num_bl, blocks, E_matrix)
lib_len = len(libraries)
print(f'Shortest path recombination complete: {lib_len} libraries found')

# Add 'M': <average chimera mutations> pairs to the libraries values.
print('Updating M', flush=True)
step2_tools.update_M(libraries, AA_alignment)
print('M updated', flush=True)

# Add 'GG_prob': <Golden Gate probability> and 'GG_sites': <valid overhangs
# for each lib_bp> pairs the libraries values.
# Requires 'normalized_ligation_counts_18h_37C.p' file to be in code directory.
print('Updating GG probabilities')
libraries = step2_tools.update_GG_prob(libraries, breakpoints)
print('GG probabilities updated')

# Convert libraries to be json-compatible and save.
libraries = [(k, v) for k, v in libraries.items()]
with open(libraries_fn, 'w') as f:
    json.dump(libraries, f)

print(f'Saved libraries to "{libraries_fn}"')
