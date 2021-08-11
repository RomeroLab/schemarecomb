import sys

from Bio import SeqIO, SeqRecord

from ggrecomb.parent_alignment import ParentAlignment
from ggrecomb.parent_alignment.pdb_structure import PDBStructure

if len(sys.argv) != 3:
    print('Usage: python trunc.py <seq_len> <num_parents>')
    sys.exit()
seq_len = int(sys.argv[1])
num_parents = int(sys.argv[2])

pa = ParentAlignment.from_fasta('../bgl3_sequences.fasta')
pdb = PDBStructure.from_pdb_file('../1GNX.pdb')
pa.pdb_structure = pdb

max_num = pa.pdb_structure.amino_acids[seq_len].res_seq_num

new_seqs = []
for i, sr in enumerate(pa.aligned_sequences[:num_parents]):
    seq = sr.seq[21:max_num].replace('-', '')
    if i == 4:
        seq = seq[2:]
    new_sr = SeqRecord.SeqRecord(seq, id=sr.id, name=sr.name,
                                 description=sr.description)
    new_seqs.append(new_sr)

pdb.amino_acids = pdb.amino_acids[:seq_len]

pa = ParentAlignment(new_seqs, pdb_structure=pdb)
print()
print('pa.align_sequences')
for sr in pa.aligned_sequences:
    print(str(sr.seq))

print()
print('pa.sequences')
for sr in pa.sequences:
    print(str(sr.seq))
print()
print('pa.pdb_structure.seq')
print(pa.pdb_structure.seq)

SeqIO.write(pa.sequences, 'trunc.fasta', 'fasta')
with open('trunc.pdb', 'w') as f:
    aas = '\n'.join([aa.to_lines() for aa in pa.pdb_structure.amino_acids])
    f.write(aas)
