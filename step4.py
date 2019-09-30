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
import sys

from Bio import SeqIO, SeqRecord, Seq

from tools import step4_tools

AA_alignment_fn = sys.argv[1]  # aligned AA sequences
CDN_alignment_fn = sys.argv[2]  # aligned codon sequences
chosen_lib_fn = sys.argv[3]  # chosen library, (lib_bps, lib_attrs)
frags_order_fn = sys.argv[4]  # fragment order output

vector_overhangs = ['TATG', 'TGAG']

# load amino acid alignment
alignment_SRs = list(SeqIO.parse(AA_alignment_fn, 'fasta'))
AA_names = [str(i.id) for i in alignment_SRs]
AA_seqs = [str(i.seq) for i in alignment_SRs]
AA_alignment = list(zip(*AA_seqs))  # list of AAs at each position

# load codon sequences
DNA_SRs = list(SeqIO.parse(CDN_alignment_fn, 'fasta'))
DNA_names = [str(i.id) for i in DNA_SRs]
assert DNA_names == AA_names
CDN_seqs = []
for sr in DNA_SRs:
    dna_seq = str(sr.seq)
    cdns = [dna_seq[index: index+3] for index in range(0, len(dna_seq), 3)]
    name = str(sr.id)
    CDN_seqs.append((name, cdns))

# load library breakpoints and overhangs
with open(chosen_lib_fn, 'r') as f:
    breakpoints, lib_attrs = json.load(f)
overhangs = lib_attrs['GG_sites']

# list of fragment to order
frags = []

# bsaI recognition sites with throwaway bases front and back
bsaI_start = 'ctagc' + 'ggtctc' + 'c'
bsaI_end = 'c' + 'gagacc' + 'gactc'

# first and last breakpoints are special -> different procedures for first and
# last blocks

# starting block
for p_name, cdn_seq in CDN_seqs:
    bp2 = breakpoints[1]  # index of breakpoint at end of block
    # make fragment sequence with breakpoint residues stripped off
    str_seq = ''.join(cdn_seq[:bp2-1]).replace('-', '')

    # sequence that gets added to front of fragment
    overhang_seq = step4_tools.first_overhang_seq(overhangs[0], str_seq)
    front_addition = bsaI_start + overhang_seq

    # sequence that gets added to back of fragment
    bp2_AA1, bp2_AA2 = AA_alignment[bp2-1], AA_alignment[bp2]
    oh2_CDN_seq = step4_tools.overhang_CDN_seq(overhangs[1], bp2_AA1, bp2_AA2,
                                               block_front_back='back')
    back_addition = oh2_CDN_seq + bsaI_end

    frag_seq = front_addition + str_seq + back_addition
    frag_sr = SeqRecord(Seq(frag_seq), p_name + '_frag1')
    frags.append(frag_sr)

# middle blocks
# breakpoint indices and overhang information for start and end of blocks
front_bps, front_ohs = breakpoints[1:-2], overhangs[1:-2]  # starts of blocks
back_bps, back_ohs = breakpoints[2:-1], overhangs[2:-1]  # ends of blocks
blocks_iter = zip(front_bps, front_ohs, back_bps, back_ohs)
for block_num, (bp1, oh1, bp2, oh2) in enumerate(2, blocks_iter):
    for p_name, cdn_seq in CDN_seqs:
        # make fragment sequence with breakpoint residues stripped off
        str_seq = ''.join(cdn_seq[bp1+1:bp2-1]).replace('-', '')

        # sequence that gets added to front of fragment
        bp1_AA1, bp1_AA2 = AA_alignment[bp1-1], AA_alignment[bp1]
        oh1_CDN_seq = step4_tools.overhang_CDN_seq(oh1, bp1_AA1, bp1_AA1,
                                                   block_front_back='front')
        front_addition = bsaI_start + oh1_CDN_seq

        # sequence that gets added to back of fragment
        bp2_AA1, bp2_AA2 = AA_alignment[bp2-1], AA_alignment[bp2]
        oh2_CDN_seq = step4_tools.overhang_CDN_seq(oh2, bp2_AA1, bp2_AA1,
                                                   block_front_back='back')
        front_addition = oh2_CDN_seq + bsaI_end

        frag_seq = front_addition + str_seq + back_addition
        frag_sr = SeqRecord(Seq(frag_seq), p_name + f'_frag{block_num}')
        frags.append(frag_sr)

# ending block
block_num += 1
for p_name, cdn_seq in CDN_seqs:
    bp1 = breakpoints[-2]  # index of breakpoint at beginning of block
    # make fragment sequence with breakpoint residues stripped off
    str_seq = ''.join(cdn_seq[bp1+1]).replace('-', '')

    # sequence that gets added to front of fragment
    bp1_AA1, bp1_AA2 = AA_alignment[bp1-1], AA_alignment[bp1]
    oh1_CDN_seq = step4_tools.overhang_CDN_seq(overhangs[-2], bp1_AA1, bp1_AA2,
                                               block_front_back='front')
    front_addition = bsaI_start + oh1_CDN_seq

    overhang_seq = step4_tools.last_overhang_seq(overhangs[-1], str_seq)
    back_addition = overhang_seq + bsaI_end

    # sequence that gets added to back of fragment
    frag_seq = front_addition + str_seq + back_addition
    frag_sr = SeqRecord(Seq(frag_seq), p_name + f'_frag{block_num}')
    frags.append(frag_sr)

SeqIO.write(frags, frags_order_fn, 'fasta')

print('Verifying fragment order against amino acid sequences', flush=True)
step4_tools.verify_fragments(frags, alignment_SRs,  breakpoints)
