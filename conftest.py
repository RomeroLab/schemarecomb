from io import StringIO
import os
from pathlib import Path
from typing import Union
from unittest import mock
from urllib import request

from Bio import SeqIO
import pytest

import ggrecomb


@pytest.fixture
def fixture_dir():
    return Path(os.path.dirname(__file__)) / 'tests' / 'fixtures/'


@pytest.fixture
def bgl3_fasta_filename(fixture_dir):
    return fixture_dir / 'bgl3_full' / 'bgl3_sequences.fasta'


@pytest.fixture
def bgl3_records(bgl3_fasta_filename):
    return list(SeqIO.parse(bgl3_fasta_filename, 'fasta'))


@pytest.fixture
def bgl3_pdb_filename(fixture_dir):
    return fixture_dir / 'bgl3_full' / '1GNX.pdb'


# Fixtures for mocking HTTP calls.


def _wrap_urlopen(**kwargs):
    """Factory for patching urlopen."""
    muscle_url = 'https://www.ebi.ac.uk/Tools/services/rest/muscle/'
    pdb_url = 'https://files.rcsb.org/'
    blast_url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

    def muscle_urlopen(muscle_request: Union[str, request.Request]):
        """Mocks MUSCLE queries on EMBI-EBI Web Services."""
        parents_jobid = 'muscle-R20210823-181055-0384-46911599-p2m'
        single_jobid = 'muscle-R20210830-210754-0277-86271462-p1m'
        if isinstance(muscle_request, request.Request):
            url = muscle_request.full_url
        else:
            url = muscle_request

        assert muscle_url in url

        if muscle_url + 'run' in url:
            # Initial run command, starts alignment job.
            from urllib.parse import parse_qs
            fasta_str, = parse_qs(muscle_request.data.decode())['sequence']
            fasta_io = StringIO(fasta_str)
            names = set(sr.name for sr in SeqIO.parse(fasta_io, 'fasta'))
            single_names = {'WP_086733846.1', 'WP_206209604.1', 'p1',
                            'WP_184902498.1', 'WP_206317985.1',
                            'WP_016827727.1'}
            parent_names = set(f'p{i}' for i in range(1, 7))
            if names == single_names:
                content = single_jobid.encode()
            elif names == parent_names:  # not single
                content = parents_jobid.encode()
            else:
                raise ValueError('Fasta query not in database.')
        elif muscle_url + 'status/' in url:
            # Status checking.
            content = b'FINISHED'
        elif muscle_url + 'result/' in url:
            # Job finished, get alignment
            assert url.split('/')[-1] == 'aln-fasta'
            jobid = url.split('/')[-2]
            if jobid == parents_jobid:
                content = kwargs['parents_aln_str'].encode()
            elif jobid == single_jobid:
                content = kwargs['single_aln_str'].encode()
            else:
                raise ValueError(f'jobid not handled: {jobid}.')

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


@pytest.fixture
def wrap_urlopen():
    return _wrap_urlopen


def _wrap_efetch(responses, acc_to_database):
    """Patches Entrez.efetch."""
    def efetch_proxy(db, id, rettype):
        assert db == 'protein'
        assert rettype == 'fasta'
        q_acc = id.split(',')[0]
        database = acc_to_database[q_acc]
        return responses[database]['fasta']
    return efetch_proxy


@pytest.fixture
def wrap_efetch():
    return _wrap_efetch


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
def bgl3_blast_SeqRecords(http_dir):
    date = database_query_dates['refseq']
    fn = http_dir / f'blast_refseq_fasta_{date}.txt'
    return list(SeqIO.parse(fn, 'fasta'))[:10]


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
def bgl3_records_aln(fixture_dir):
    fn = fixture_dir / 'bgl3_full' / 'bgl3_sequences_aln.fasta'
    return list(SeqIO.parse(fn, 'fasta'))


@pytest.fixture
def bgl3_parents_aln_str(bgl3_records_aln):
    seqs_f = StringIO('')
    SeqIO.write(bgl3_records_aln, seqs_f, 'fasta')
    seqs_f.seek(0)
    return seqs_f.read()


@pytest.fixture
def bgl3_single_aln_str(fixture_dir):
    fn = fixture_dir / 'bgl3_single/' / 'found_parents_aln.fasta'
    records = list(SeqIO.parse(fn, 'fasta'))
    seqs_f = StringIO('')
    SeqIO.write(records, seqs_f, 'fasta')
    seqs_f.seek(0)
    return seqs_f.read()


@pytest.fixture
def mock_bgl3_blast_query(bgl3_records, mocker, blast_http_responses,
                          rid_to_database, acc_to_database, wrap_urlopen,
                          wrap_efetch, bgl3_parents_aln_str,
                          bgl3_single_aln_str, bgl3_pdb_filename):
    """Wrapper over wrap_urlopen and wrap_efetch with preexisting args.

    This:
        fake_urlopen, fake_efetch = mock_bgl3_blast_query
        mocker.patch('ggrecomb.parent_alignment.urlopen', fake_urlopen)
        mocker.patch('ggrecomb.parent_alignment.Entrez.efetch', fake_efetch)

    is equivalent to this:
        query_seq = str(bgl3_records[0].seq)
        refseq_responses = blast_http_responses('refseq')
        pdb_responses = blast_http_responses('pdb')
        responses = {'refseq_protein': refseq_responses,
                     'pdbaa': pdb_responses}
        mocker.patch(
            'ggrecomb.parent_alignment.urlopen',
            wrap_urlopen(
                query_seq=query_seq,
                parents_aln_str=bgl3_parents_aln_str,
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

    """
    query_seq = str(bgl3_records[0].seq)
    refseq_responses = blast_http_responses('refseq')
    pdb_responses = blast_http_responses('pdb')
    responses = {'refseq_protein': refseq_responses, 'pdbaa': pdb_responses}
    fake_urlopen = wrap_urlopen(
               query_seq=query_seq,
               parents_aln_str=bgl3_parents_aln_str,
               single_aln_str=bgl3_single_aln_str,
               responses=responses,
               rid_to_database=rid_to_database,
               bgl3_pdb_filename=bgl3_pdb_filename,
    )
    fake_efetch = wrap_efetch(
            responses=responses,
            acc_to_database=acc_to_database
        )
    return fake_urlopen, fake_efetch


@pytest.fixture
def bgl3_mock_namespace(doctest_namespace, mocker, mock_bgl3_blast_query):
    fake_urlopen, fake_efetch = mock_bgl3_blast_query
    mocker.patch('ggrecomb.parent_alignment.urlopen', fake_urlopen)
    mocker.patch('ggrecomb.parent_alignment.Entrez.efetch', fake_efetch)
