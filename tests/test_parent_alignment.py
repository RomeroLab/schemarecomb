from itertools import combinations

from Bio import SeqIO
import pytest

from ggsr import parent_alignment
from ggsr.parent_alignment.blast_query import blast_query
from ggsr.parent_alignment.utils import iden_diff


def mock_urlopen(request):
    url = request._full_url
    data = dict(x.split('=') for x in request._data.decode().split('&'))
    if url == 'https://blast.ncbi.nlm.nih.gov/Blast.cgi':
        # initial BLAST query
        if all([data.get('DATABASE', None) == 'refseq_protein',
                data['CMD'] == 'Put']):
            return open('bgl3_sample/blast_put_20210622.txt', 'rb')

        # check if complete, return completion right away (not true behavior)
        elif data['CMD'] == 'Get' and 'RESULTS_FILE' not in data:
            return open('bgl3_sample/blast_status-20_20210622.txt', 'rb')

        # get accession numbers from BLAST request
        elif data['CMD'] == 'Get' and 'RESULTS_FILE' in data:
            return open('bgl3_sample/blast_get_20210622.txt', 'rb')

    raise ValueError('Handling for this request is not implemented.')


def mock_efetch(db, id, rettype):
    if db == 'protein' and rettype == 'fasta':  # id is way too long
        return open('bgl3_sample/fasta_str_20210622.txt')
    raise ValueError('Handling for this request is not implemented.')


ParentAlignment = parent_alignment.ParentAlignment
parent_alignment.blast_query.urlopen = mock_urlopen
parent_alignment.blast_query.Entrez.efetch = mock_efetch


def check_aligned_sequences(pa, intended_sequences_len):
    assert len(pa.sequences) == intended_sequences_len
    assert len(pa.aligned_sequences) == intended_sequences_len

    for sr, aln_sr in zip(pa.sequences, pa.aligned_sequences):
        assert sr.seq == aln_sr.seq.replace('-', '')

    seq1 = pa.aligned_sequences[0].seq
    assert all(len(seq1) == len(aln_sr.seq) for aln_sr in pa.aligned_sequences)


def compare_seqrecords(sr1, sr2):
    assert sr1.id == sr2.id
    assert sr1.name == sr2.name
    assert sr1.description == sr2.description
    assert str(sr1.seq) == str(sr2.seq)


def test_no_input():
    with pytest.raises(ValueError):
        ParentAlignment([])


def test_from_fasta():  # tests obtain_seqs and choose_candidates also
    pa = ParentAlignment.from_fasta('bgl3_sample/bgl3_sequences.fasta')
    original_seqs = pa.sequences
    check_aligned_sequences(pa, 6)
    pa.obtain_seqs(7)
    check_aligned_sequences(pa, 7)
    for sr1, sr2 in zip(original_seqs, pa.sequences[:6]):
        compare_seqrecords(sr1, sr2)


def test_from_fasta_no_align():
    pa = ParentAlignment.from_fasta('bgl3_sample/bgl3_sequences.fasta',
                                    auto_align=False)
    original_seqs = pa.sequences
    assert pa.aligned_sequences is None
    pa.obtain_seqs(7)
    assert pa.aligned_sequences is None
    pa.align()
    check_aligned_sequences(pa, 7)
    for sr1, sr2 in zip(original_seqs, pa.sequences[:6]):
        compare_seqrecords(sr1, sr2)


def test_from_single():
    in_sr = list(SeqIO.parse('bgl3_sample/bgl3_sequences.fasta', 'fasta'))[0]
    ParentAlignment.from_single(in_sr, num_final_sequences=6,
                                desired_identity=0.65)
    in_string = str(in_sr.seq)
    in_name = in_sr.id
    pa = ParentAlignment.from_single(in_string, in_name)
    assert pa.sequences[0].id == pa.sequences[0].name

    with pytest.raises(ValueError):
        ParentAlignment.from_single(in_string)


def test_json():
    p_aln = ParentAlignment.from_fasta('bgl3_sample/bgl3_sequences.fasta')
    aln_json = p_aln.to_json()
    p_aln2 = ParentAlignment.from_json(aln_json)
    vars1 = vars(p_aln)
    vars2 = vars(p_aln2)
    assert vars1.keys() == vars2.keys()
    for k, v in vars(p_aln).items():
        if k in ('sequences', 'aligned_sequences'):
            for sr1, sr2 in zip(v, vars2[k]):
                compare_seqrecords(sr1, sr2)
        else:
            assert v == vars2[k]


def combo_choose_cands(sorted_cand_diffs, target_identity, num_choose):
    best_diff = 1.0
    best_sets = []
    for srs in combinations(sorted_cand_diffs, num_choose):
        max_diff = max([sr[1] for sr in srs])
        if max_diff > best_diff:
            continue
        combos = combinations([sr[0] for sr in srs], 2)
        max_cc = max(iden_diff(sr1, sr2, target_identity) for sr1, sr2
                     in combos)
        max_diff = max(max_diff, max_cc)
        if max_diff < best_diff:
            best_diff = max_diff
            best_sets = [set(sr.id for sr, _ in srs)]
        elif max_diff == best_diff:
            best_sets.append(set(sr.id for sr, _ in srs))
    return best_diff, best_sets


def test_choose_candidates():
    bgl3_seqs = list(SeqIO.parse('bgl3_sample/bgl3_sequences.fasta',
                                 'fasta'))
    num_start = 3
    num_choose = 2
    num_final = num_start + num_choose
    in_seqs = bgl3_seqs[:num_start]
    pa = ParentAlignment(in_seqs, False)

    # Test 1: choose from small set, auto-calc indentity.
    cands = bgl3_seqs[num_start:]
    identity = pa._check_identity(None)
    pa.add_from_candidates(cands, num_final_sequences=num_final)

    cand_diffs = []
    for cand in cands:
        max_diff = max(iden_diff(cand, x, identity) for x in pa.sequences)
        cand_diffs.append((cand, max_diff))
    sorted_cand_diffs = list(sorted(cand_diffs, key=lambda x: x[1]))
    _, combo_sets = combo_choose_cands(sorted_cand_diffs, identity, num_choose)

    pa_chosen = pa.sequences[3:]
    assert set(pa.id for pa in pa_chosen) in combo_sets

    # Test 2: ask for too many seqs
    pa.add_from_candidates([], 1)

    # Test 3: choose from BLAST query, 0.6 identity.
    query_sr = pa.sequences[0]
    query_seq = str(query_sr.seq)
    cands = list(blast_query(query_seq))

    cand_diffs = []
    for cand in cands:
        max_diff = max(iden_diff(cand, x, 0.6) for x in pa.sequences)
        cand_diffs.append((cand, max_diff))

    sorted_cand_diffs = list(sorted(cand_diffs, key=lambda x: x[1]))
    _, combo_sets = combo_choose_cands(sorted_cand_diffs, 0.6, 2)

    pa.add_from_candidates(cands, 7, 0.6)
    pa_chosen = pa.sequences[num_final:]  # num_final left over from Test 1

    assert set(pa.id for pa in pa_chosen) in combo_sets
