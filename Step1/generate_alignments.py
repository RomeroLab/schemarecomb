import os
from subprocess import Popen

from Bio import Seq, SeqRecord, SeqIO

import tools

in_fn = 'sequences.fasta'        # name of sequences file read by muscle
out_file_prefix = 'alignment'    # prefix of output files
pdb_fn = '1GNX.pdb'

# Align sequences with muscle
Popen(f'muscle -in {in_fn} -out {out_file_prefix}_AA.fasta',
      shell=True).wait()

alignment_SRs = list(SeqIO.parse(f'{out_file_prefix}_AA.fasta', 'fasta'))
names = [str(i.id) for i in alignment_SRs]
descriptions = [str(i.description) for i in alignment_SRs]
AA_seqs = [str(i.seq) for i in alignment_SRs]
AA_alignment = list(zip(*AA_seqs))  # get list of AAs at each position

pdb_AAs = tools.read_PDB(pdb_fn)
pdb_Seq = Seq.Seq(tools.get_PDB_seq(pdb_AAs))
pdb_SeqRecord = SeqRecord.SeqRecord(pdb_Seq, name='PDB',
                                    description=f'{pdb_fn} sequence')
temp_SRs = [pdb_SeqRecord] + alignment_SRs
SeqIO.write(temp_SRs, 'temp_seqs.fasta', 'fasta')
Popen(f'muscle -in temp_seqs.fasta -out PDB_align.fasta', shell=True).wait()
os.delete('temp_seqs.fasta')
pdb_align_SRs = SeqIO.read('PDB_align.fasta', 'fasta')
pdb_align_seqs = [str(i.seq) for i in pdb_align_SRs]
pdb_alignment = list(zip(*AA_seqs))

pdb_residues = iter(pdb_AAs)
pos = 0
for pdb_AA, parent_AAs in pdb_alignment:
    valid_parent_pos = all([a == '-' for a in parent_AAs])
    if valid_parent_pos:
        pos += 1
    if pdb_AA != '-':
        AA = next(pdb_residues)
        assert tools.seq1(AA.resName) == pdb_AA
        if valid_parent_pos:
            AA.resSeq = pos
        else:
            AA.resSeq = None

# Add stop codons if necessary
if AA_alignment[-1] != tuple(['*'] * len(alignment_SRs)):
    AA_alignment.append(tuple(['*'] * len(alignment_SRs)))


CDN_align = [tools.get_pos_CDNs(pos) for pos in AA_alignment]
CDN_seqs = [''.join(CDN_pos) for CDN_pos in zip(*CDN_align)]
# save_SeqRecords(CDN_seqs, f'{out_file_prefix}_CDN.fasta')
