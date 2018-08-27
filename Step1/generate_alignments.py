from subprocess import Popen
from Bio import SeqIO
import pickle
import random

AA_C31 = {'A': ('GCT', 'GCA'), 'R': ('CGT', 'CGA'), 'N': ('AAT',), 'D': ('GAT',), 'C': ('TGT',), 'Q': ('CAA', 'CAG'), 'E': ('GAA',), 'G': ('GGT',), 'H': ('CAT', 'CAC'), 'I': ('ATT', 'ATC'), 'L': ('TTA', 'TTG', 'CTA'), 'K': ('AAA',), 'M': ('ATG',), 'F': ('TTT',), 'P': ('CCT', 'CCA'), 'S': ('AGT', 'TCA'), 'T': ('ACA', 'ACT'), 'W': ('TGG',), 'Y': ('TAT',), 'V': ('GTT', 'GTA'), '*': ('TGA',), '-': ('---',)}

in_fn = 'sequences.fasta'        #name of sequences read by muscle
out_file_prefix = 'alignment'             #prefix of output files

Popen('muscle -in {} -out {}_AA.fasta'.format(in_fn, out_file_prefix), shell=True).wait()
alignment_file = list(SeqIO.parse('{}_AA.fasta'.format(out_file_prefix), 'fasta'))

pickle.dump([str(i.id) for i in alignment_file], open('{}_names.p'.format(out_file_prefix), 'wb'))

AA_alignment = list(zip(*[str(i.seq) for i in alignment_file]))
AA_alignment.append(tuple(['*'] * len(alignment_file)))            #add stop AAs to each sequence
pickle.dump(AA_alignment, open('{}_AA.p'.format(out_file_prefix), 'wb'))

single_AAs = [set(pos) for pos in AA_alignment]
single_CDNs = [{AA: random.choice(AA_C31[AA]) for AA in pos} for pos in single_AAs]
CDN_align = [tuple(CDNs[AA] for AA in align_pos) for align_pos, CDNs in zip(AA_alignment, single_CDNs)]
pickle.dump(CDN_align, open('{}_CDN.p'.format(out_file_prefix), 'wb'))
