import time
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.request import Request

from Bio import Entrez, SeqIO
from Bio.Blast.NCBIWWW import _parse_qblast_ref_page
# ^This import is Python heresy but who's going to stop me?

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def blast_query(query_seq, maximum_hits=10000):
    '''Gets sequences from query_seq BLAST hits.

    Direct copy of Peter Cock's Bio.Blast.NCBIWWW.qblast, with simplifications
    and added progress print-outs.

    Returns:
        Fasta generator for sequence hits.
    '''
    params = [
        ('PROGRAM', 'blastp'),
        ('DATABASE', 'refseq_protein'),
        ('QUERY', query_seq),
        ('HITLIST_SIZE', maximum_hits),
        ("CMD", "Put"),
    ]

    message = urlencode(params).encode()
    request = Request(NCBI_BLAST_URL, message, {"User-Agent":
                                                "BiopythonClient"})
    handle = urlopen(request)

    rid, rtoe = _parse_qblast_ref_page(handle)

    print(rid)

    params = [
        ('RID', rid),
        ('FORMAT_TYPE', 'XML'),
        ('CMD', 'Get'),
    ]
    message = urlencode(params).encode()

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

    print(results)

    print('Job complete.')

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

    headers = handle.readline().decode().strip().split(',')
    acclen_index = headers.index('Acc. Len')
    acc_index = -1  # split by ',' unfortunately splits last field in two
    candidate_accessions = []
    for line in handle:
        fields = line.decode().strip().split(',')
        if len(fields) <= 1:
            continue
        accession = fields[acc_index].replace('"', '').replace(')', '')
        candidate = fields[acclen_index], accession
        candidate_accessions.append(candidate)

    # TODO: filter based on acclen
    # TODO: query for new sequences on demand (generator)
    # TODO: Error handling and input verification throughout module

    accessions = [x[1] for x in candidate_accessions]
    acc_str = ','.join(accessions)
    Entrez.email = 'bbremer@bu.edu'  # this probably needs to change
    fasta_str = Entrez.efetch(db='protein', id=acc_str, rettype='fasta')
    
    print('returning')

    return SeqIO.parse(fasta_str, 'fasta')
