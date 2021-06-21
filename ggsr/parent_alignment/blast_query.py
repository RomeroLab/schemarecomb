# blast_query.py

"""Wrapper for querying BLAST over the web."""

import time
from collections.abc import Generator
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.request import Request

from Bio import Entrez, SeqIO, SeqRecord
from Bio.Blast.NCBIWWW import _parse_qblast_ref_page
# ^This import is Python heresy but who's going to stop me?

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def blast_query(query_seq: str, database: str = 'refseq_protein',
                maximum_hits: int = 10000) -> Generator[SeqRecord.SeqRecord]:
    """Gets sequences from query_seq BLAST hits.

    Direct copy of Peter Cock's Bio.Blast.NCBIWWW.qblast, with simplifications
    and added progress print-outs.

    Args:
        query_seq: Query sequence for BLAST search.
        database: BLAST database to search. Must be either 'refseq_protein' or
            'pdbaa'
        maximum_hits: Maximum number of hits allowed to be returned. Affects
            runtime: reducing this number is faster but you risk missing good
            candidate sequences.

    Returns:
        FASTA generator for sequence hits.
    """

    # Only refseq_protein and pdbaa databases supported.
    if database not in ('refseq_protein', 'pdbaa'):
        raise ValueError('database must be "refseq_protein" or "pdbaa"')

    # Initiate BLAST search.
    params = [
        ('PROGRAM', 'blastp'),
        ('DATABASE', database),
        ('QUERY', query_seq),
        ('HITLIST_SIZE', maximum_hits),
        ("CMD", "Put"),
    ]
    message = urlencode(params).encode()
    request = Request(NCBI_BLAST_URL, message, {"User-Agent":
                                                "BiopythonClient"})
    handle = urlopen(request)

    # Get BLAST Request ID.
    rid, _ = _parse_qblast_ref_page(handle)

    # Setup BLAST job checking request.
    params = [
        ('RID', rid),
        ('FORMAT_TYPE', 'XML'),
        ('CMD', 'Get'),
    ]
    message = urlencode(params).encode()

    # Query BLAST once a minute until job is done.
    previous = 0
    delay = 20  # seconds
    while True:
        current = time.time()
        wait = previous + delay - current
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current
        # delay by at least 60 seconds only if running the request against the
        # public NCBI API
        if delay < 60:
            # Wasn't a quick return, must wait at least a minute
            delay = 60

        print('Checking for job completion...')
        request = Request(NCBI_BLAST_URL, message, {"User-Agent":
                                                    "BiopythonClient"})
        handle = urlopen(request)
        results = handle.read().decode()

        # Can see an "\n\n" page while results are in progress,
        # if so just wait a bit longer...
        if results == "\n\n":
            continue
        # XML results don't have the Status tag when finished
        if "Status=" not in results:
            break
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i + len("Status="): j].strip()
        if status.upper() == "READY":
            break

    print('Job complete.')

    # Request to get hit accessions from search results.
    params = [
        ('RID', rid),
        ('RESULTS_FILE', 'on'),
        ('FORMAT_TYPE', 'CSV'),
        ('DESCRIPTIONS', maximum_hits),
        ('FORMAT_OBJECT', 'Alignment'),
        ('DOWNLOAD_TEMPL', 'Results'),
        ('QUERY_INDEX', '0'),
        ('CMD', 'Get'),
    ]
    message = urlencode(params).encode()
    request = Request(NCBI_BLAST_URL, message, {"User-Agent":
                                                "BiopythonClient"})
    handle = urlopen(request)

    # Unpack accessions into candidate_accessions list.
    handle.readline()  # headers
    candidate_accessions = []
    for line in handle:
        fields = line.decode().strip().split(',')
        if len(fields) <= 1:
            continue
        # Split by ',' unfortunately splits last field in two.
        accession = fields[-1].replace('"', '').replace(')', '')
        candidate_accessions.append(accession)

    # TODO: query for new sequences on demand (generator)

    # Get accession sequences from NCBI.
    acc_str = ','.join(candidate_accessions)
    Entrez.email = 'bbremer@bu.edu'  # this probably needs to change
    fasta_str = Entrez.efetch(db='protein', id=acc_str, rettype='fasta')

    print('returning')

    return SeqIO.parse(fasta_str, 'fasta')
