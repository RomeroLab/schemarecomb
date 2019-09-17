import os

from Bio import Seq, SeqRecord, SeqIO

import tools

in_fn = 'test_seqs.fasta'        # name of sequences file read by muscle
out_file_prefix = 'alignment'    # prefix of output files
pdb_fn = 'test.pdb'

tools.run_muscle(in_fn, out_file_prefix + '_AA.fasta')

# Load muscle input fasta file
alignment_SRs = list(SeqIO.parse(f'{out_file_prefix}_AA.fasta', 'fasta'))
names = [str(i.id) for i in alignment_SRs]
descriptions = [str(i.description) for i in alignment_SRs]
AA_seqs = [str(i.seq) for i in alignment_SRs]
AA_alignment = list(zip(*AA_seqs))  # get list of AAs at each position

# Load sequence from pdb and temporarily add to AA alignment
pdb_AAs = tools.read_PDB(pdb_fn)
pdb_Seq = tools.get_PDB_seq(pdb_AAs)
pdb_Seq = Seq.Seq(pdb_Seq)
pdb_SeqRecord = SeqRecord.SeqRecord(pdb_Seq, name='PDB',
                                    description=f'{pdb_fn} sequence')
temp_SRs = [pdb_SeqRecord] + alignment_SRs

# Align pdb sequence to AA aligned sequences
SeqIO.write(temp_SRs, 'temp_seqs.fasta', 'fasta')
tools.run_muscle('temp_seqs.fasta', 'PDB_align.fasta')
os.remove('temp_seqs.fasta')
pdb_align_SRs = SeqIO.parse('PDB_align.fasta', 'fasta')
pdb_align_seqs = [str(i.seq) for i in pdb_align_SRs]
pdb_alignment = list(zip(*pdb_align_seqs))

# Renumber pdb residues and save new pdb
pdb_AAs_residue = iter(pdb_AAs)
pos = 0    # PDB will be 1-indexed for consistency with other PDB files
for pdb_AA, *parent_AAs in pdb_alignment:
    valid_parent_pos = any([a != '-' for a in parent_AAs])
    if valid_parent_pos:
        pos += 1
    if pdb_AA != '-':
        AA = next(pdb_AAs_residue)
        assert tools.seq1(AA.resName) == pdb_AA
        if valid_parent_pos:
            AA.resSeq = pos
        else:
            AA.resSeq = None
tools.write_PDB(pdb_AAs, 'renumbered_' + pdb_fn)

# Calculate contacts
pdb_AAs = list(tools.read_PDB('renumbered_' + pdb_fn))
for AA in pdb_AAs:
    AA.resSeq -= 1   # make index start from 0
contacts = []
for index, AA1 in enumerate(pdb_AAs):
    for AA2 in pdb_AAs[index+1:]:
        if AA1.d(AA2) < 4.5:
            contacts.append((AA1.resSeq, AA2.resSeq))

# Add stop codons if necessary
if AA_alignment[-1] != tuple(['*'] * len(alignment_SRs)):
    AA_alignment.append(tuple(['*'] * len(alignment_SRs)))

# Make and save codon alignments
CDN_align = [tools.get_pos_CDNs(pos) for pos in AA_alignment]
CDN_seqs = [''.join(CDN_pos) for CDN_pos in zip(*CDN_align)]
tools.save_SeqRecords(CDN_seqs, names, descriptions, 
                      f'{out_file_prefix}_CDN.fasta')
