import copy
from io import StringIO
from typing import Union
from unittest import mock
from urllib import request

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pytest

import ggrecomb
from ggrecomb import ParentSequences
from ggrecomb import PDBStructure

ggrecomb  # Just to get rid of import not used warning.


@pytest.fixture
def bgl3_records_aln(fixture_dir):
    fn = fixture_dir / 'bgl3_full' / 'bgl3_sequences_aln.fasta'
    return list(SeqIO.parse(fn, 'fasta'))


@pytest.fixture
def bgl3_aln_str(bgl3_records_aln):
    seqs_f = StringIO('')
    SeqIO.write(bgl3_records_aln, seqs_f, 'fasta')
    seqs_f.seek(0)
    return seqs_f.read()


@pytest.fixture
def bgl3_PDBStructure(bgl3_pdb_filename):
    return PDBStructure.from_pdb_file(bgl3_pdb_filename)


def equal(obj1: ParentSequences, obj2: ParentSequences) -> bool:
    if obj1.__class__ is not obj2.__class__:
        print('class')
        return False

    if obj1.__dict__.keys() != obj2.__dict__.keys():
        print('keys')
        return False

    for k in obj1.__dict__:
        if k == 'records':
            recs1 = obj1.records
            recs2 = obj2.records
            for r1, r2 in zip(recs1, recs2):
                if str(r1.seq) != str(r2.seq) or r1.name != r2.name \
                   or r1.id != r2.id or r1.description != r2.description:
                    print('records')
                    return False
        elif k == 'pdb_structure':
            if obj1.pdb_structure.seq != obj2.pdb_structure.seq:
                print('pdb_structure')
                return False
        else:
            if not obj1.__dict__[k] == obj2.__dict__[k]:
                print('key:', k)
                return False

    return True


def wrap_urlopen(**kwargs):
    """Factory for patching urlopen."""
    muscle_url = 'https://www.ebi.ac.uk/Tools/services/rest/muscle/'
    pdb_url = 'https://files.rcsb.org/'
    blast_url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

    def muscle_urlopen(url: Union[str, request.Request]):
        """Mocks MUSCLE queries on EMBI-EBI Web Services."""
        jobid = 'muscle-R20210823-181055-0384-46911599-p2m'

        if isinstance(url, request.Request):
            url = url.full_url

        assert muscle_url in url

        if muscle_url + 'run' in url:
            # Initial run command, starts alignment job.
            content = jobid.encode()
        elif url == muscle_url + 'status/' + jobid:
            # Status checking.
            content = b'FINISHED'
        elif url == muscle_url + 'result/' + jobid + '/aln-fasta':
            # Job finished, get alignment.
            content = kwargs['aln_str'].encode()

        def fake_read():
            return content

        http_response = mock.Mock()
        http_response.read = fake_read

        return http_response

    def pdb_urlopen(url: Union[str, request.Request]):
        """Mocks Protein Data Bank queries."""
        # TODO: Finish this
        if isinstance(url, str):
            if url == 'https://files.rcsb.org/view/1GNX.pdb':
                return open(kwargs['bgl3_pdb_filename'], 'rb')
        else:
            url = url.full_url
        raise ValueError('PDB url not handled: ' + url)

    def blast_urlopen(url: Union[str, request.Request]):
        """Mocks BLAST webtool."""
        if isinstance(url, str):
            raise ValueError('blast_urlopen input must be a Request')

        req_data = url._data.decode()
        data_dict = {}
        for item in req_data.split('&'):
            k, v = item.split('=')
            data_dict[k] = v
        url = url.full_url

        database_query_size = {'refseq_protein': '10000', 'pdbaa': '100'}

        if data_dict['CMD'] == 'Put':
            # BLAST init
            # Only these options supported in mock.
            assert data_dict['PROGRAM'] == 'blastp'
            assert data_dict['QUERY'] == kwargs['query_seq']

            database = data_dict['DATABASE']
            assert data_dict['HITLIST_SIZE'] == database_query_size[database]

            return kwargs['responses'][database]['put']

        if data_dict['CMD'] == 'Get' and data_dict['FORMAT_TYPE'] == 'XML':
            # Status.
            database = kwargs['rid_to_database'][data_dict['RID']]
            return kwargs['responses'][database]['status']

        if data_dict['CMD'] == 'Get' and data_dict['FORMAT_TYPE'] == 'CSV':
            # Get accessions.
            database = kwargs['rid_to_database'][data_dict['RID']]
            assert data_dict['RESULTS_FILE'] == 'on'
            assert data_dict['DESCRIPTIONS'] == database_query_size[database]
            assert data_dict['FORMAT_OBJECT'] == 'Alignment'
            assert data_dict['DOWNLOAD_TEMPL'] == 'Results'
            assert data_dict['QUERY_INDEX'] == '0'
            return kwargs['responses'][database]['accessions']

        # Not handled
        raise ValueError('Unhandled blast query passed in.')

    def urlopen_proxy(url: Union[str, request.Request]):
        """Based on url, call the appropriate mock function."""
        # Get url as string.
        temp_url = url
        if isinstance(temp_url, request.Request):
            temp_url = temp_url.full_url

        # Return the right url.
        if muscle_url in temp_url:
            return muscle_urlopen(url)
        elif pdb_url in temp_url:
            return pdb_urlopen(url)
        elif blast_url in temp_url:
            return blast_urlopen(url)
        else:
            raise ValueError('URL not handled: ' + temp_url)

    return urlopen_proxy


def wrap_efetch(responses, acc_to_database):
    """Patches Entrez.efetch."""
    def efetch_proxy(db, id, rettype):
        assert db == 'protein'
        assert rettype == 'fasta'
        q_acc = id.split(',')[0]
        database = acc_to_database[q_acc]
        return responses[database]['fasta']
    return efetch_proxy


@pytest.fixture
def http_dir(fixture_dir):
    return fixture_dir / 'bgl3_full' / 'http'


# For easy updating if queries are remade.
database_query_dates = {'refseq': '20210622', 'pdb': '20210628'}


@pytest.fixture
def rid_to_database(http_dir):
    """Maps BLAST RID to database from that run using initial http response."""
    r_t_d = {}
    for database, date in database_query_dates.items():
        fn = http_dir / f'blast_{database}_put_{date}.txt'
        with open(fn, 'rb') as handle:
            rid, _ = ggrecomb.parent_alignment._parse_qblast_ref_page(handle)
        if database == 'refseq':
            # Patch because refseq filename does not match database query name
            database = 'refseq_protein'
        elif database == 'pdb':
            database = 'pdbaa'
        r_t_d[rid] = database
    return r_t_d


@pytest.fixture
def acc_to_database(http_dir):
    """Map first sequence accession from blast accessions files to database."""
    a_t_d = {}
    for database, date in database_query_dates.items():
        fn = http_dir / f'blast_{database}_accessions_{date}.txt'
        with open(fn, 'rb') as handle:
            handle.readline()
            first_acc_line = handle.readline()
            fields = first_acc_line.decode().strip().split(',')
            acc = fields[-1].replace('"', '').replace(')', '')
        if database == 'refseq':
            # Patch because refseq filename does not match database query name
            database = 'refseq_protein'
        elif database == 'pdb':
            database = 'pdbaa'
        a_t_d[acc] = database
    return a_t_d


@pytest.fixture
def blast_http_responses(http_dir):
    """Sim HTTP responses from files in http_dir."""
    def get_response(database, request_type):
        """Given database and request_type, yield mock HTTP response file."""
        date = database_query_dates[database]
        fn = http_dir / f'blast_{database}_{request_type}_{date}.txt'
        if request_type == 'fasta':
            yield open(fn)  # Entrez return in text mode.
        else:
            yield open(fn, 'rb')

    def get_database_responses(database):
        """Get the four blast HTTP response type files in dict."""
        f_dict = {}
        for request_t in 'put', 'status', 'accessions', 'fasta':
            f_dict[request_t] = next(get_response(database, request_t))
        return f_dict

    return get_database_responses


@pytest.fixture
def bgl3_blast_SeqRecords(http_dir):
    date = database_query_dates['refseq']
    fn = http_dir / f'blast_refseq_fasta_{date}.txt'
    return list(SeqIO.parse(fn, 'fasta'))[:10]


def test_query_blast(bgl3_records, mocker, blast_http_responses,
                     rid_to_database, acc_to_database):
    query_seq = str(bgl3_records[0].seq)

    # patch urlopen for speed
    refseq_responses = blast_http_responses('refseq')
    pdb_responses = blast_http_responses('pdb')
    responses = {'refseq_protein': refseq_responses, 'pdbaa': pdb_responses}
    mocker.patch(
        'ggrecomb.parent_alignment.urlopen',
        wrap_urlopen(
            query_seq=query_seq,
            aln_str=bgl3_aln_str,
            responses=responses,
            rid_to_database=rid_to_database
        )
    )
    # also need to patch Entrez.efetch for BLAST runs
    mocker.patch(
        'ggrecomb.parent_alignment.Entrez.efetch',
        wrap_efetch(
            responses=responses,
            acc_to_database=acc_to_database
        )
    )

    blast_fastas = list(ggrecomb.parent_alignment.query_blast(query_seq))
    assert blast_fastas
    assert all(isinstance(bf, SeqRecord) for bf in blast_fastas)

    blast_pdb_fastas = list(
        ggrecomb.parent_alignment.query_blast(query_seq, 'pdbaa', 100)
    )
    assert blast_pdb_fastas
    assert all(isinstance(bpf, SeqRecord) for bpf in blast_pdb_fastas)


def test_choose_candidates(bgl3_records, bgl3_blast_SeqRecords):
    # TODO: This is a simple test. Make it more robust for 0.2.0.
    preexisting = bgl3_records[:2]
    desired_identity = 0.7
    num_additional = 2

    # Brute force way
    from itertools import combinations
    from itertools import product
    calc_identity = ggrecomb.parent_alignment.calc_identity
    combos = combinations(bgl3_blast_SeqRecords, num_additional)
    diff_cands = []
    for candidates in combos:
        max_diff = max(abs(desired_identity - calc_identity(cand, preex))
                       for cand, preex in product(candidates, preexisting))
        cand_diff = abs(desired_identity - calc_identity(*candidates))
        max_diff = max(max_diff, cand_diff)
        cand_names = set(cand.name for cand in candidates)
        diff_cands.append((max_diff, cand_names))
    # Collect all optimal sets.
    opt_diff = min(dc[0] for dc in diff_cands)
    brute_cands = [names for diff, names in diff_cands if diff == opt_diff]

    # Using choose_candidates.
    choose_cands = ggrecomb.parent_alignment.choose_candidates(
        bgl3_blast_SeqRecords,
        preexisting,
        num_additional,
        desired_identity
    )

    # choose_candidates answer should be in brute force optimals.
    assert set(c.name for c in choose_cands) in brute_cands


class TestParentSequences:
    def test_init_invalid_input(self, bgl3_records):
        # No sequences.
        with pytest.raises(ValueError):
            ParentSequences([])

        # Conflicting auto_alignment, prealignment.
        with pytest.raises(ValueError):
            ParentSequences(bgl3_records, auto_align=True, prealigned=True)

    def test_align(self, bgl3_records, mocker, bgl3_aln_str, bgl3_records_aln):
        mocker.patch('ggrecomb.parent_alignment.urlopen',
                     wrap_urlopen(aln_str=bgl3_aln_str))
        parents = ParentSequences(bgl3_records)
        parents.align()

        # Test against alignment we expect to see.
        exp_alignment = zip(*[str(rec.seq) for rec in bgl3_records_aln])
        assert all(p == a for p, a in zip(parents.alignment, exp_alignment))

        # Test resulting ParentAlignment against prealigned.
        parents_aln = ParentSequences(bgl3_records_aln, prealigned=True)
        assert parents.alignment == parents_aln.alignment
        # But the records are different.
        assert parents.p0_aligned != parents_aln.alignment

    @pytest.mark.parametrize(
        'pdb,auto_align,prealigned',
        [
            (False, True, False),
            (False, False, True),
            (False, False, False),
            (True, True, False),
            (True, False, True),
            (True, False, False),
            # (*, True, True) in test_init_invalid_input
        ]
    )
    def test_init(self, bgl3_records, bgl3_records_aln, bgl3_PDBStructure,
                  pdb, auto_align, prealigned, mocker, bgl3_aln_str):
        if pdb:
            input_pdb = bgl3_PDBStructure
        else:
            input_pdb = None

        if prealigned:
            input_records = bgl3_records_aln
        else:
            input_records = bgl3_records

        if auto_align:
            # align method will be called, patch the HTTP call.
            mocker.patch('ggrecomb.parent_alignment.urlopen',
                         wrap_urlopen(aln_str=bgl3_aln_str))

        parents = ParentSequences(input_records, input_pdb, auto_align,
                                  prealigned)

        assert parents.records
        for rec, prec in zip(input_records, parents.records):
            assert rec is prec

        # If either one is true, parents should have an alignment after
        # initialization.
        hasaln = prealigned or auto_align
        if hasaln:
            assert parents.alignment
            aligned_seqs = [''.join(seq) for seq
                            in zip(*parents.alignment)]
            # Make sure aligned_seqs matches records.
            for rec, seq in zip(parents.records, aligned_seqs):
                assert str(rec.seq).replace('-', '') == seq.replace('-',
                                                                    '')

        else:
            # No alignment set during init.
            with pytest.raises(AttributeError):
                parents.alignment
            with pytest.raises(AttributeError):
                parents.p0_aligned

        if input_pdb is None:
            # No pdb_structure set during init.
            with pytest.raises(AttributeError):
                parents.pdb_structure
        else:
            pdb_struct = parents.pdb_structure
            if hasaln:
                # Make sure pdb_structure was renumbered and matches
                # alignment.
                assert parents.p0_aligned == pdb_struct.renumbering_seq
                assert pdb_struct.seq.replace('-', '')
                assert all(p1_aa == pdb_aa for p1_aa, pdb_aa
                           in zip(parents.p0_aligned, pdb_struct.seq)
                           if pdb_aa != '-')
            else:
                with pytest.raises(AttributeError):
                    # pdb_struct is not renumbered.
                    pdb_struct.renumbering_seq

    def test_attr_pdb(self, bgl3_records, bgl3_records_aln, bgl3_PDBStructure,
                      mocker, bgl3_aln_str):
        """Test attribute changes starting with pdb."""
        # Unaligned, no PDB.
        parents = ParentSequences(bgl3_records)
        assert parents.records
        with pytest.raises(AttributeError):
            parents.alignment
        with pytest.raises(AttributeError):
            parents.pdb_structure

        # Add pdb_structure, should not be renumbered.
        parents.pdb_structure = bgl3_PDBStructure
        with pytest.raises(AttributeError):
            parents.alignment
        assert parents.pdb_structure
        with pytest.raises(AttributeError):
            self.pdb_structure.renumbering_seq

        # Direct alignment is not allowed.
        with pytest.raises(AttributeError):
            parents.alignment = 'dummy string'
        with pytest.raises(AttributeError):
            parents.alignment

        aln_seqs = [str(sr.seq) for sr in bgl3_records_aln]

        # Calling align method, patch HTTP call.
        mocker.patch('ggrecomb.parent_alignment.urlopen',
                     wrap_urlopen(aln_str=bgl3_aln_str))

        # align(), then new_alignment.
        parents.align()
        assert parents.p0_aligned == parents.pdb_structure.renumbering_seq
        parents_copy = copy.deepcopy(parents)
        parents.new_alignment(aln_seqs)
        assert equal(parents, parents_copy)

        parents.new_alignment(aln_seqs)
        assert parents.p0_aligned == parents.pdb_structure.renumbering_seq
        parents_copy = copy.deepcopy(parents)
        parents.align()
        assert equal(parents, parents_copy)

    def test_attr_alignment(self, bgl3_records_aln, bgl3_PDBStructure):
        """Test attribute changes starting with alignment."""
        # Aligned, no PDB.
        parents = ParentSequences(bgl3_records_aln, prealigned=True)
        parents_unchanged = copy.deepcopy(parents)

        assert parents.alignment
        with pytest.raises(AttributeError):
            parents.pdb_structure

        previous_alignment = copy.deepcopy(parents.alignment)
        parents.pdb_structure = bgl3_PDBStructure
        assert previous_alignment == parents.alignment
        assert parents.pdb_structure
        assert parents.pdb_structure.renumbering_seq  # pdb is aligned
        assert parents.p0_aligned == parents.pdb_structure.renumbering_seq

        del parents.pdb_structure
        assert equal(parents, parents_unchanged)

    def test_add_from_candidates(self, bgl3_records, bgl3_blast_SeqRecords):
        parents = ParentSequences(bgl3_records)

        with pytest.raises(ValueError):
            parents.add_from_candidates([], 7, 0.7)
        with pytest.raises(ValueError):
            parents.add_from_candidates(bgl3_blast_SeqRecords, 7, 0.0)
        with pytest.raises(ValueError):
            parents.add_from_candidates(bgl3_blast_SeqRecords, 7, 1.0)

        parents.add_from_candidates(bgl3_blast_SeqRecords, 8, 0.7)
        assert len(parents.records) == 8
        assert all(isinstance(sr, SeqRecord) for sr in parents.records)

    def test_obtain_seqs(self, bgl3_records, mocker, bgl3_blast_SeqRecords):
        # We already tested list(query_blast(query_seq)) with the bgl3
        # query_seq in test_query_blast, so we can patch it to instantly return
        # the result.
        mocker.patch('ggrecomb.parent_alignment.query_blast',
                     return_value=bgl3_blast_SeqRecords)

        parents = ParentSequences(copy.deepcopy(bgl3_records))
        parents.obtain_seqs(8, 0.7)

        # default identity
        parents = ParentSequences(bgl3_records)
        parents.obtain_seqs(8)

    def test_get_PDB(self, bgl3_records, mocker, fixture_dir,
                     blast_http_responses, rid_to_database, acc_to_database):
        # We already tested list(query_blast(query_str, 'pdbaa', 100)) with
        # the bgl3 query_seq in test_query_blast, so we can patch it to
        # instantly return the result.
        query_seq = str(bgl3_records[0].seq)

        # patch urlopen for speed
        refseq_responses = blast_http_responses('refseq')
        pdb_responses = blast_http_responses('pdb')
        responses = {
            'refseq_protein': refseq_responses,
            'pdbaa': pdb_responses
        }
        bgl3_pdb_filename = fixture_dir / 'bgl3_full' / '1GNX.pdb'
        mocker.patch(
            'ggrecomb.parent_alignment.urlopen',
            wrap_urlopen(
                query_seq=query_seq,
                aln_str=bgl3_aln_str,
                responses=responses,
                rid_to_database=rid_to_database,
                bgl3_pdb_filename=bgl3_pdb_filename
            )
        )
        # also need to patch Entrez.efetch for BLAST runs
        mocker.patch(
            'ggrecomb.parent_alignment.Entrez.efetch',
            wrap_efetch(
                responses=responses,
                acc_to_database=acc_to_database
            )
        )

        parents = ParentSequences(bgl3_records)
        with pytest.raises(AttributeError):
            parents.pdb_structure
        parents.get_PDB()
        assert isinstance(parents.pdb_structure, PDBStructure)

    def test_json(self, bgl3_records, bgl3_PDBStructure, mocker, bgl3_aln_str):

        # mock urlopen for muscle call
        mocker.patch(
            'ggrecomb.parent_alignment.urlopen',
            wrap_urlopen(
                aln_str=bgl3_aln_str,
            )
        )

        # No pdb_structure or alignment.
        in_parents = ParentSequences(bgl3_records)
        json_str = in_parents.to_json()
        out_parents = ParentSequences.from_json(json_str)
        assert equal(in_parents, out_parents)

        # Alignment but no pdb_structure.
        in_parents = ParentSequences(bgl3_records, auto_align=True)
        json_str = in_parents.to_json()
        out_parents = ParentSequences.from_json(json_str)
        assert equal(in_parents, out_parents)

        # pdb_structure but no alignment.
        in_parents = ParentSequences(bgl3_records, bgl3_PDBStructure)
        json_str = in_parents.to_json()
        out_parents = ParentSequences.from_json(json_str)
        assert equal(in_parents, out_parents)

        # Alignment but no pdb_structure.
        in_parents = ParentSequences(bgl3_records, bgl3_PDBStructure,
                                     auto_align=True)
        json_str = in_parents.to_json()
        out_parents = ParentSequences.from_json(json_str)
        assert equal(in_parents, out_parents)
