import random

from Bio import Seq, SeqRecord, SeqIO
from Bio.SeqUtils import seq1

# List of optimal codons for E. coli.
AA_C31 = {'A': ('GCT', 'GCA'), 'R': ('CGT', 'CGA'), 'N': ('AAT',),
          'D': ('GAT',), 'C': ('TGT',), 'Q': ('CAA', 'CAG'), 'E': ('GAA',),
          'G': ('GGT',), 'H': ('CAT', 'CAC'), 'I': ('ATT', 'ATC'),
          'L': ('TTA', 'TTG', 'CTA'), 'K': ('AAA',), 'M': ('ATG',),
          'F': ('TTT',), 'P': ('CCT', 'CCA'), 'S': ('AGT', 'TCA'),
          'T': ('ACA', 'ACT'), 'W': ('TGG',), 'Y': ('TAT',),
          'V': ('GTT', 'GTA'), '*': ('TGA',), '-': ('---',)}


class atom(object):
    def __init__(self, line):
        self.resName = line[17:20].strip()
        self.resSeq = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())


def d(a1, a2):
    return ((a1.x - a2.x)**2 + (a1.y - a2.y)**2 + (a1.z - a2.z)**2) ** (1/2)


class residue(object):
    def __init__(self, atom):
        self.atoms = [atom]
        self.resName = atom.resName
        self._resSeq = atom.resSeq

    def add_atom(self, atom):
        self.atoms.append(atom)

    def d(self, res):
        m = d(self.atoms[0], res.atoms[0])
        for i in self.atoms:
            for j in res.atoms:
                dis = d(i, j)
                if dis < m:
                    m = dis
        return m

    @property
    def resSeq(self):
        return self._resSeq

    @resSeq.setter
    def resSeq(self, resSeq):
        for a in self.atoms:
            a.resSeq = resSeq
        self.resSeq = resSeq


def read_PDB(pdb_fn):
    with open(pdb_fn) as f:
        atoms = [atom(i.strip()) for i in f if 'ATOM' in i[:4]]
    AAs = {}
    for a in atoms:
        if a.resSeq not in AAs:
            AAs.update({a.resSeq: residue(a)})
        else:
            AAs[a.resSeq].add_atom(a)
    return AAs.values()


def get_PDB_seq(AAs):
    return [seq1(AA.resName) for AA in AAs]


def save_SeqRecords(seqs, names, descriptions, handle):
    seqs = [Seq.Seq(sequence) for sequence in seqs]
    zipped_SRs = zip(seqs, names, descriptions)
    records = [SeqRecord.SeqRecord(s, name=n, description=d)
               for s, n, d in zipped_SRs]
    SeqIO.write(records, handle, 'fasta')


def get_pos_CDNs(pos):
    single_AAs = set(pos)
    translation = {AA: random.choice(AA_C31[AA]) for AA in single_AAs}
    return tuple(translation[AA] for AA in pos)
