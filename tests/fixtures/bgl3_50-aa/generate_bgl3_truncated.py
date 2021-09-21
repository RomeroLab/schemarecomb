"""Truncates bgl3_full fixtures for faster runtime.

The runtime of this script is slow (use web services) because it is only run
once.
"""

import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from schemarecomb import ParentSequences
from schemarecomb import PDBStructure

trunc_dir = os.path.dirname(os.path.realpath(__file__))
bgl3_full_dir = os.path.join(trunc_dir, '..', 'bgl3_full')

seq_len = 50  # final size of alignment between pdb and p0_aligned
num_parents = 6

pdb_fn = os.path.join(bgl3_full_dir, '1GNX.pdb')
pdb = PDBStructure.from_pdb_file(pdb_fn)
parents_fn = os.path.join(bgl3_full_dir, 'bgl3_sequences.fasta')
parents = ParentSequences.from_fasta(parents_fn, pdb_structure=pdb)

# Reduce number of parents.
parents.records = parents.records[:num_parents]

# Get parent alignment and renumber pdb.
print(f'aligning {num_parents} parents')
parents.align()

# Find the first index that doesn't have any gaps in the alignment.
for pdb_index, aa in enumerate(parents.pdb_structure.amino_acids):
    aln_index = aa.index
    if '-' not in parents.alignment[aln_index]:
        break
first_aa_index = aln_index
last_aa_index = parents.pdb_structure.amino_acids[pdb_index+seq_len-1].index

print(first_aa_index, last_aa_index)

new_records = []
aligned_sequences = [''.join(seq) for seq in zip(*parents.alignment)]
for sr, aln_seq in zip(parents.records, aligned_sequences):
    new_seq = aln_seq[first_aa_index:last_aa_index+1]
    new_sr = SeqRecord(Seq(new_seq), id=sr.id, name=sr.name,
                       description=sr.description)
    new_records.append(new_sr)

reduced_parents = ParentSequences(new_records, pdb, prealigned=True)

# Write parents
new_parents_fn = os.path.join(trunc_dir, 'bgl3_trunc_aln.fasta')
SeqIO.write(reduced_parents.records, new_parents_fn, 'fasta')

# Write pdb
with open('1GNX_trunc.pdb', 'w') as f:
    pdb = reduced_parents.pdb_structure
    aa_strs = ['\n'.join(aa.to_lines()) for aa in pdb.amino_acids]
    out_str = '\n'.join(aa_strs)
    f.write(out_str)
