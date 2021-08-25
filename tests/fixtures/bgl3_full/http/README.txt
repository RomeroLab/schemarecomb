Data from previous HTTP for testing purposes.

Filenames follow this conventions:
<url ident>_<database>_<request type>_<access date>.txt

blast put files are the blast initiation HTTP request:
>>> params = [
>>>    ('PROGRAM', 'blastp'),
>>>    ('DATABASE', database),
>>>    ('QUERY', query_seq),
>>>    ('HITLIST_SIZE', maximum_hits),
>>>    ("CMD", "Put"),
>>> ]
>>> message = urlencode(params).encode()
>>> Request("https://blast.ncbi.nlm.nih.gov/Blast.cgi", message)
where query_seq is the 1GNX_T10M sequence.

blast status files are the blast status check:
>>> params = [
>>>     ('RID', rid),
>>>     ('FORMAT_TYPE', 'XML'),
>>>     ('CMD', 'Get'),
>>> ]
>>> message = urlencode(params).encode()
>>> Request("https://blast.ncbi.nlm.nih.gov/Blast.cgi", message)
Where rid is found from the response to the initiation HTTP request.
Note that blast status files have the "finished" response for fast testing
purposes, BLAST also returns "in progress" responses to this query if the
BLAST run is not done.

blast accessions files obtain the hit sequence names:
>>> params = [
>>>     ('RID', rid),
>>>     ('RESULTS_FILE', 'on'),
>>>     ('FORMAT_TYPE', 'CSV'),
>>>     ('DESCRIPTIONS', maximum_hits),
>>>     ('FORMAT_OBJECT', 'Alignment'),
>>>     ('DOWNLOAD_TEMPL', 'Results'),
>>>     ('QUERY_INDEX', '0'),
>>>     ('CMD', 'Get'),
>>> ]
>>> message = urlencode(params).encode()
>>> request = Request(NCBI_BLAST_URL, message, {"User-Agent":
...                                             "BiopythonClient"})
>>> handle = urlopen(request)

blast fasta files are the hit FASTAs return from Entrez:
>>> acc_str = ','.join(candidate_accessions)
>>> Entrez.email = '<user email>'  # this probably needs to change
>>> fasta_str = Entrez.efetch(db='protein', id=acc_str, rettype='fasta')
where acc_str is the comma-separated hit accessions from blast accessions.
