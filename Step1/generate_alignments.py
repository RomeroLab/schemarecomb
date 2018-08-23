from subprocess import Popen
from Bio import SeqIO
import pickle
import random

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T','W', 'Y', 'V','*', '-']

C31 = [['GCT','GCA'],['CGT','CGA'],['AAT'],['GAT'],['TGT'],['CAA','CAG'],['GAA'],['GGT'],['CAT','CAC'],['ATT','ATC'],['TTA','TTG','CTA'],['AAA'],['ATG'],['TTT'],['CCT','CCA'],['AGT','TCA'],['ACA','ACT'],['TGG'],['TAT'],['GTT','GTA'],['TGA'],['---']]

in_fn = 'sequences'
out_fn = 'alignment'

Popen('muscle -in {}.fasta -out {}.fasta'.format(in_fn, out_fn), shell=True).wait()
alignment_file = list(SeqIO.parse('{}.fasta'.format(out_fn), 'fasta'))

pickle.dump([str(i.id) for i in alignment_file], open('names_{}.p'.format(out_fn), 'wb'))

AA_alignment = list(zip(*[str(i.seq) for i in alignment_file]))
AA_alignment.append(tuple(['*'] * len(alignment_file)))
pickle.dump(AA_alignment, open('AA_{}.p'.format(out_fn), 'wb'))

single_AAs = [set(i) for i in AA_alignment]
single_CDNs = [{j:random.choice(C31[AA.index(j)]) for j in i} for i in single_AAs]
CDN_align = [tuple([single_CDNs[ii][j] for j in i]) for ii,i in enumerate(AA_alignment)]
pickle.dump(CDN_align, open('CDN_{}.p'.format(out_fn), 'wb'))
