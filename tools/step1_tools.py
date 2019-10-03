import random
from subprocess import Popen
import sys
import urllib

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


def get_pdb(sequence):
    """ Find and download the closestPDB sequence to sequence.

    Uses a BLAST search of the PDB database to find the PDB sequence closest to
    the input sequence. The BLAST search sequenceIdentityCutoff is scanned from
    0 to 100 with a binary search algorithm until a unique pdb is found. If
    no unique PDB is found, but one exists, the returned PDBs all have the same
    percent identity with sequence and any can be chosen.

    Documentation for PDB API: https://www.rcsb.org/pages/webservices. It
    would have been more efficient to use the direct BLAST-PDB API
    (https://www.rcsb.org/pdb/software/rest.do#blastPdb) to return BLAST search
    sequence identities, but it did not return the best sequences for a bgl3
    example. I may have been doing something wrong or that API is broken.

    Args:
        sequence: Query sequence for BLAST search

    Returns:
        pdb filename of (possibly one of) the database sequence with highest
            pairwise identity to sequence, where pdb filename was downloaded
    """

    url = 'http://www.rcsb.org/pdb/rest/search'
    lower_bound, upper_bound = 0, 100  # adjust cutoff with average of these
    while True:
        identity_cutoff = (upper_bound + lower_bound) // 2
        print("Current identity percentage: %s" % identity_cutoff, flush=True)

        # pdb-blast search parameters
        queryText = """
        <orgPdbQuery>
            <queryType>org.pdb.query.simple.SequenceQuery</queryType>
            <sequence>%s</sequence>
            <searchTool>blast</searchTool>
            <maskLowComplexity>yes</maskLowComplexity>
            <sequenceIdentityCutoff>%s</sequenceIdentityCutoff>
        </orgPdbQuery>""" % (sequence, identity_cutoff)

        # make request and read output
        req = urllib.request.Request(url, data=queryText.encode())
        f = urllib.request.urlopen(req)
        result = f.read().decode()
        pdbs = result.strip().split('\n')

        # if unique pdb found, select it and break
        if (len(pdbs) == 1 and pdbs != ['']):
            pdb = pdbs[0]
            break

        # if upper_bound and lower_bound are close, break (2 since 1 can lead
        # to infinite loop)
        if upper_bound - lower_bound < 2:
            break

        # need another refinement round, if no pdbs, avg is too high. if more
        # than 1 pdb, avg is too low, save a pdb for possible break next round
        if pdbs == ['']:
            upper_bound = identity_cutoff
        elif len(pdbs) > 1:
            lower_bound = identity_cutoff
            pdb = pdbs[0]

    if not pdb:
        print('no pdb found')
        sys.exit()
    else:
        print(f'pdb found: {pdb}')

    # pdb is format "1GNX:1", not sure what the ":1" is but remove it
    pdb = pdb.split(':')[0]

    # retrieve the pdb file
    pdb_fn = pdb + '.pdb'
    url = 'http://www.rcsb.org/pdb/files/' + pdb_fn
    with urllib.request.urlopen(url) as response, \
            open(pdb_fn, 'wb') as out_file:
        data = response.read()
        out_file.write(data)

    return pdb_fn


class atom(object):
    def __init__(self, line):
        self.line = line
        self.resName = line[17:20].strip()
        self.resSeq = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())

    def __str__(self):
        line = list(self.line)
        line[22:26] = list('{:4d}'.format(self.resSeq))
        return ''.join(line)


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
    def resSeq(self, new_resSeq):
        for a in self.atoms:
            a.resSeq = new_resSeq
        self._resSeq = new_resSeq

    def __repr__(self):
        return f'residue: {self.resName} {self.resSeq}'


def run_muscle(in_fn, out_fn):
    # Align sequences with muscle
    print(f'Running muscle with {in_fn}, saving to {out_fn}.', flush=True)
    Popen(f'muscle -in {in_fn} -out {out_fn}',
          shell=True).wait()
    out_SRs = {SR.name: SR for SR in SeqIO.parse(out_fn, 'fasta')}
    in_labels = [SR.name for SR in SeqIO.parse(in_fn, 'fasta')]
    stable_out_SRs = [out_SRs[label] for label in in_labels]
    SeqIO.write(stable_out_SRs, out_fn, 'fasta')
    print('\n\nMuscle done.')


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


def write_PDB(residues, fn):
    with open(fn, 'w') as fn:
        for AA in residues:
            if AA.resSeq is not None:
                for atom in AA.atoms:
                    fn.write(str(atom) + '\n')


def get_PDB_seq(AAs):
    return ''.join([seq1(AA.resName) for AA in AAs])


def save_SeqRecords(seqs, names, handle):
    seqs = [Seq.Seq(sequence) for sequence in seqs]
    zipped_SRs = zip(seqs, names)
    records = [SeqRecord.SeqRecord(s, id=n) for s, n in zipped_SRs]
    SeqIO.write(records, handle, 'fasta')


def get_pos_CDNs(pos):
    single_AAs = set(pos)
    translation = {AA: random.choice(AA_C31[AA]) for AA in single_AAs}
    return tuple(translation[AA] for AA in pos)
