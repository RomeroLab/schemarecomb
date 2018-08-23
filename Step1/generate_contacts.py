from Bio.SeqUtils import seq1
from Bio import SeqIO
import pickle

pdb_fn = '5ng6_a.pdb'
pdb_seq_name_in_alignment = '5ng6_a'
alignment_fn = 'alignment'

class atom:
    def __init__(self, line):
        self.serial = line[6:11].strip()
        self.name = line[12:16].strip()
        self.altLoc = line[16].strip()
        self.resName = line[17:20].strip()
        self.chainID = line[21].strip()
        self.resSeq = int(line[22:26].strip())
        self.iCode = line[26].strip()
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.element = line[76:78]
        self.charge = line[78:80]
  
def d(a1, a2):
    return ((a1.x - a2.x)**2 + (a1.y - a2.y)**2 + (a1.z - a2.z)**2) ** (1/2)    
    
class residue:
    def __init__(self, atom):
        self.atoms = [atom]
        self.resName = atom.resName
        
    def add_atom(self,atom):
        self.atoms.append(atom)
        
    def d(self, res):
        m = d(self.atoms[0], res.atoms[0])
        for i in self.atoms:
            for j in res.atoms:
                dis = d(i, j)
                if dis < m:
                    m = dis
        return m

with open(pdb_fn) as f:
    lines = [i.strip() for i in f.readlines() if 'ATOM' in i[:4]]

atoms = [atom(line) for line in lines]
AAs = {}
for a in atoms:
    if a.resSeq not in AAs:
        AAs.update({a.resSeq:residue(a)})
    else:
        AAs[a.resSeq].add_atom(a)

initial_contacts = []
for i in sorted(AAs):
    print('\rCurrent AA: {}'.format(i), end=' ', flush=True)
    for j in sorted(AAs)[sorted(AAs).index(i)+1:]:
        if AAs[i].d(AAs[j]) < 4.5:
            initial_contacts.append([i,j])
print()

pdb_index_in_alignment = pickle.load(open('names_{}.p'.format(alignment_fn), 'rb')).index(pdb_seq_name_in_alignment)
wt_msa = str(list(SeqIO.parse('AA_{}.fasta'.format(alignment_fn), 'fasta'))[pdb_index_in_alignment].seq)
wt_pdb_ex = ''.join([seq1(AAs[i+1].resName) if i+1 in AAs else '-' for i in range(len(wt_msa))])

PDB_index = 0
MSA_index = 0
while True:
    print('\nCurrent PDB_index: {}, Current MSA_index: {}'.format(PDB_index, MSA_index))
    print('PDB seq: ' + wt_pdb_ex[PDB_index : PDB_index + 50] + '...')
    print('MSA seq: ' + wt_msa[MSA_index : MSA_index + 50] + '...')
    i = input('Are these sequences aligned? [Y/N]: ')
    while i != 'Y' and i != 'N':
        i = input('Please enter Y or N: ')
    if i == 'Y':
        break
    else:
        PDB_index = int(input('Enter new PDB_ind: '))
        MSA_index = int(input('Enter new MSA_ind: '))

old_new = {}
while PDB_index < len(wt_pdb_ex) and MSA_index < len(wt_msa):
    if wt_msa[MSA_index] != '-' and wt_pdb_ex[PDB_index] != '-':   
        old_new.update({PDB_index+1:MSA_index})   #convert back to 1-indexing of PDB file
        PDB_index += 1
        MSA_index += 1
    elif wt_pdb_ex[PDB_index] == '-' and wt_msa[MSA_index] != '-':
        PDB_index += 1
        MSA_index += 1
    elif wt_msa[MSA_index] == '-':
        MSA_index += 1

for k,v in old_new.items():
    assert seq1(AAs[k].resName) == wt_msa[v]

new_contacts = {(old_new[i[0]], old_new[i[1]]):1 for i in initial_contacts}

pickle.dump(new_contacts, open('contacts.p', 'wb'))
