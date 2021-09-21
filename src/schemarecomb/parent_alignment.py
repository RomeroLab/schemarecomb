# parent_alignment.py

"""Alignment of parental sequences for protein recombination.

This module provides the definition of :class:`schemarecomb.ParentSequences`,
which represents an alignment of parental protein sequences used in
schemarecomb library creation. The module can be used directly to use accessory
functions for calculating identity between sequences, querying the BLAST web
interface, or selecting additional parent sequences from a list of candidates.

"""

from collections import defaultdict
from itertools import combinations
from io import StringIO
import json
import os
from subprocess import CalledProcessError
from subprocess import run
import time
from typing import Optional
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib.request import Request

from Bio import Entrez
from Bio import pairwise2
from Bio import SeqIO
from Bio.Blast.NCBIWWW import _parse_qblast_ref_page
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import schemarecomb


def calc_identity(sr1: SeqRecord, sr2: SeqRecord) \
        -> float:
    """Calculate the BLAST identity between two sequences.

    Note:
        "BLAST identity" is defined as the number of matches in a pairwise
        alignment divided by the length of alignment. This is different from
        the "traditional" definition of identity: number of matches divided by
        the minimum of the sequence lengths.

    Parameters:
        sr1: SeqRecord of first sequence to compare.
        sr2: SeqRecord of second sequence to compare.

    Returns:
        The BLAST identity between the two sequences.

    """
    s1 = str(sr1.seq)
    s2 = str(sr2.seq)
    aln_score = pairwise2.align.globalxx(s1, s2, score_only=True)
    aln_len = len(s1) + len(s2) - aln_score
    return aln_score / aln_len


def query_blast(
    query_seq: str,
    database: str = 'refseq_protein',
    maximum_hits: int = 10000
) -> SeqIO.FastaIO.FastaIterator:
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

    Raises:
        ValueError: If database is not "refseq_protein" or "pdbaa".
    """

    NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

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
    # TODO: Change from BiopythonClient.
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
    request = Request(NCBI_BLAST_URL, message, {"User-Agent":
                                                "BiopythonClient"})

    # Query BLAST once a minute until job is done.
    previous = 0.0
    delay = 20.0  # seconds
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

        # print('Checking for job completion...')
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

    # print('Job complete.')

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
    Entrez.email = 'bbremer@wisc.edu'  # this probably needs to change
    fasta_str = Entrez.efetch(db='protein', id=acc_str, rettype='fasta')

    # print('returning')

    return SeqIO.parse(fasta_str, 'fasta')


class _NoNewNodeError(RuntimeError):
    """Raised when no node is added to the tree."""
    pass


class _TreeNode:
    """Represents a set of candidate sequences.

    Helper class for choose_candidates function.

    Attributes:
        cands: List of candidates SeqRecords.
        max_cc_diff: Maximum out of abs(% identity - target_identity) between
        each pair of candidates.
        max_pc_diff: Maximum out of abs(% identity - target_identity) of each
            parent and each candidate.
        children: Child nodes of this node.
        parent: Parent node of this node.
    """
    def __init__(
        self,
        cands: SeqRecord = [],
        max_cc_diff: float = 0.0,
        max_pc_diff: float = 0.0,
        parent: '_TreeNode' = None
    ) -> None:
        """Initialize node instance with precomputed identies."""
        self.cands = cands
        self.max_cc_diff = max_cc_diff
        self.max_pc_diff = max_pc_diff
        self.children: list['_TreeNode'] = []
        if parent is not None:
            self.parent = parent

    def new_cand(self, cand: SeqRecord, new_pc_diff: float,
                 diff_thresh: float, target_identity: float) -> '_TreeNode':
        """Create new child node with candidates self.cands + [cand].

        Args:
            cand: New candidate sequence.
            new_pc_diff: Maximum out of abs(% identity - target_identity)
                between cand and each parent. This is assumed to be greater
                than the "pc_diff" of anything in self.cands, by virtue of
                pre-sorting based on the pc_diff.
            diff_thresh: Current known optimum for a valid set of candidates.
                If the potentially new set has a higher diff than this value,
                we can skip creation because this node and all its children
                are non-optimal.
            target_identity: Ideal cross-wise identity between all sequences in
                the concatenation set of parents and selected candidates.
        """
        if self.cands:
            # Calculate the identity diff between the new candidate and each
            # candidate in self.cands.
            new_cc_diff = max(abs(target_identity - calc_identity(cand, x))
                              for x in self.cands)
            # Choose maximum between old cc_diff and new_cc_diff.
            new_cc_diff = max(self.max_cc_diff, new_cc_diff)
        else:
            new_cc_diff = 0.0  # only one candidate in new node
        if max(new_cc_diff, new_pc_diff) >= diff_thresh:
            # don't create a new node because it is already non-optimal
            raise _NoNewNodeError
            # return None  # probably change b/c EP2 Item 20
        new_cands = self.cands + [cand]
        return _TreeNode(new_cands, new_cc_diff, new_pc_diff, self)

    def add_cand(self, cand: SeqRecord, new_pc_diff: float,
                 diff_thresh: float, target_identity: float) -> None:
        """new_cand new_cand for adding non-leaf nodes."""
        try:
            new_node = self.new_cand(cand, new_pc_diff, diff_thresh,
                                     target_identity)
        except _NoNewNodeError:
            return
        if new_node is not None:
            self.children.append(new_node)

    def delete(self) -> None:
        """Delete node from tree."""
        self.parent.children.remove(self)

    def __repr__(self) -> str:
        """Return CSV string of the node's candidate ids."""
        return ','.join([c.id for c in self.cands])


class _Tree:
    """Data structure for finding the ideal set of candidate sequences.

    Helper class for choose_candidates function.

    Attributes:
        base: Node representing the empty set of candidates.
        num_seq: Number of sequences to choose. Depth of tree is num_seq + 1.
        target_identity: Ideal cross-wise identity between all sequences in
            the concatenation set of parents and selected candidates.
        best_leaf: The TreeNode with the current best set of candidates.
        best_diff: The maximum abs(% identity - target_identity) out of each
            pair in best_leaf.cands + parents.
    """
    def __init__(self, num_seqs: int, target_identity: float):
        """Tree for finding <num_seqs> candidates with <target_identity>."""
        self.base = _TreeNode()
        self.num_seqs = num_seqs
        self.target_identity = target_identity
        self.best_diff = 1.0

    def add_cand(self, cand: SeqRecord, new_pc_diff: float) -> Optional[str]:
        """Add cand to candidate sets in tree and add new nodes if needed.

        Args:
            cand: New candidate sequence.
            new_pc_diff: Maximum out of abs(% identity - target_identity)
                between cand and each parent. This is assumed to be greater
                than the "pc_diff" of anything in self.cands, by virtue of
                pre-sorting based on the pc_diff.
        """
        # Breadth-first traversal over tree, adding new cand to each node.
        stack = [self.base]
        while stack:
            curr_node = stack.pop()

            # Prune tree if curr_node's max identity difference is greater
            # than the current known best.
            curr_max = max(curr_node.max_cc_diff, curr_node.max_pc_diff)
            if curr_max > self.best_diff:
                curr_node.delete()
                continue

            # Leaf node case: if new node will have num_seqs cands, then we can
            # evaluate it directly and don't need to add it to tree.
            if len(curr_node.cands) == self.num_seqs - 1:
                try:
                    leaf = curr_node.new_cand(cand, new_pc_diff,
                                              self.best_diff,
                                              self.target_identity)
                except _NoNewNodeError:
                    continue
                # leaf max diff is smaller than self.best_diff, new best!
                self.best_leaf = leaf
                self.best_diff = max(leaf.max_cc_diff, leaf.max_pc_diff)

                # Every future leaf will have equal or greater max_pc_diff,
                # so this leaf must be the best one and we can stop.
                if leaf.max_pc_diff >= leaf.max_cc_diff:
                    return 'best found'

                continue

            # Non-leaf node case.
            stack += curr_node.children  # comes first b/c don't want dup cands
            curr_node.add_cand(cand, new_pc_diff, self.best_diff,
                               self.target_identity)
        return None

    def shape(self) -> list[int]:
        """Number of nodes at each level in tree."""
        shape: dict[int, int] = defaultdict(int)
        stack = [self.base]
        while stack:
            curr_node = stack.pop()
            shape[len(curr_node.cands)] += 1
            stack += curr_node.children
        return [shape[i] for i, _ in enumerate(shape)]


def choose_candidates(
    candidate_sequences: list[SeqRecord],
    existing_parents: list[SeqRecord],
    num_additional: int,
    desired_identity: float
) -> list[SeqRecord]:
    """Choose the ideal set of candidate sequences.

    Ideal set is defined as the set of candidates with the minimum max_diff,
    where max_diff is the maximum cross-wise abs(% identity - desired_identity)
    between each pair of sequences in the concatenation of the set and parents.

    Args:
        candidate_sequences: Sequences able to be selected.
        existing_parents: Parent sequences in the library already.
        num_additional: Number of candidate sequences to choose.
        desired_identity: Ideal cross-wise identity between all sequences in
            the concatenation set of parents and selected candidates.

    Returns:
        Ideal set of candidate SeqRecords.

    Raises:
        ValueError: if nonpositive num_additional provided, desired_identity
            not bet 0.0 and 1.0 (noninclusive), or less candidate_sequences
            than num_additional provided.
    """

    # Parameter checking.
    if num_additional < 1:
        raise ValueError('num_additional must be positive.')
    if not 0.0 < desired_identity < 1.0:
        raise ValueError('desired_identity must be between 0.0 and 1.0, '
                         'exclusive.')
    if len(candidate_sequences) < num_additional:
        raise ValueError('Insufficient number of candidates provided.')
    if not existing_parents:
        raise ValueError('existing_parents must not be empty.')
    # TODO: Handle empty existing_parents as valid? I.e. choose starting from
    # nothing. (Version 0.2.0)

    # For each candidate, find the maximum identity difference with each
    # parent.

    cand_diffs = []
    # TODO: logging (Version 0.2.0)
    # print('Number of candidates:', len(candidate_sequences))
    # print('Calculating pc_diff and sorted')
    for i, cand in enumerate(candidate_sequences):
        # print(i, '\r', end='')
        max_diff = max(abs(desired_identity - calc_identity(cand, x))
                       for x in existing_parents)
        cand_diffs.append((cand, max_diff))
    # print()

    # Sort candidates by difference to enable short-circuiting. sorted is fast
    # enough, but this could be faster with heapq.
    sorted_cand_diffs = list(sorted(cand_diffs, key=lambda x: x[1]))

    # Construct Tree and find the best set. In the worst case this might take
    # awhile.
    # TODO: Bound the time this requires?
    # print('Constructing tree.')
    tree = _Tree(num_additional, desired_identity)
    for i, (cand, pc_diff) in enumerate(sorted_cand_diffs):
        # print(i, '\r', end='')
        if pc_diff > tree.best_diff:
            break
        ret = tree.add_cand(cand, pc_diff)
        if ret == 'best found':
            break
    # print()
    return tree.best_leaf.cands


def web_muscle(records: list[SeqRecord]) -> str:
    """Align sequences using MUSCLE on EMBL-EBI Web Services.

    Parameters:
        records: The SeqRecords holding the sequences to align.

    Return:
        The alignment in a FASTA format string. It can be read manually or
            with "SeqIO.parse(StringIO(aln_str), 'fasta')".

    """
    email = 'bbremer@wisc.edu'
    base_url = 'https://www.ebi.ac.uk/Tools/services/rest/muscle/'

    seqs_f = StringIO('')
    SeqIO.write(records, seqs_f, 'fasta')
    seqs_f.seek(0)
    seqs_str = seqs_f.read()

    # Submit job.
    params = {
        'email': email,
        'format': 'fasta',
        'sequence': seqs_str,
    }
    message = urlencode(params).encode()
    url = base_url + 'run'
    request = Request(url, message)
    jobid = urlopen(request).read().decode()

    # Wait for completion.
    url = base_url + 'status/' + jobid
    response = urlopen(url).read().decode()
    while response == 'RUNNING':
        time.sleep(1)
        response = urlopen(url).read().decode()
    if response != 'FINISHED':
        url = base_url + 'result/' + jobid + '/error'
        error = urlopen(url).read().decode()
        raise RuntimeError(f"MUSCLE status response: {repr(response)}, "
                           f"error message: {error}")

    # Get resulting alignment.
    url = base_url + 'result/' + jobid + '/aln-fasta'
    aln_str = urlopen(url).read().decode()
    return aln_str


def local_muscle(records: list[SeqRecord]) -> str:
    """Align sequences with MUSCLE and set to aligned_sequences.

    WARNING: Experimental in the current version.
    """
    # print('running alignment')

    IN_FN = 'temp_muscle_input.fasta'
    OUT_FN = 'temp_muscle_output.fasta'

    # make file for MUSCLE input
    # TODO: Use tempfile package.
    SeqIO.write(records, IN_FN, 'fasta')

    # run MUSCLE
    try:
        run(f'muscle -in {IN_FN} -out {OUT_FN}', shell=True, check=True)
    except CalledProcessError:
        # print('Something is wrong with MUSCLE call. Is MUSCLE installed?')
        raise
    finally:
        os.remove(IN_FN)

    with open(OUT_FN) as f:
        aln_str = f.read()
    os.remove(OUT_FN)
    return aln_str


class _ParentSequences:
    """
    Parent protein sequences for recombinant library design.

    This class sets up the data needed to run recombinant design algorithms,
    e.g. intaking and aligning parental sequences, aligning the parents, and
    finding a PDB structure or additional parents. Instances of this class are
    passed into functions further down the schemarecomb pipeline.

    Note:
        The first sequence (records[0]) has special importance, as it's used to
        find PDB structures or additional parents. For parent sequence sets
        with high enough identity (~60%), this should generally not be an
        issue, but you may get different results by changing the order of the
        sequences.

    Parameters:
        records: Parental amino acid sequences for schemarecomb
            calculations.
        auto_align: If True, the records are aligned upon initialization.
        prealigned: If True, the records are already aligned upon
            initialization. auto_align and prealigned must not both be
            True.
        pdb_structure: Protein Data Bank structure that represents the three-
            dimensional structure of the aligned parent sequences. Will be
            automatically renumbered if either auto_align or prealigned is
            True.

    Attributes:
        records (list[SeqRecord]): Parental amino acid sequences for
            schemarecomb calculations. `BioPython SeqRecords
            <https://biopython.org/wiki/SeqRecord>`_ contain sequence metadata
            such as the name of the sequence. The (unaligned) ith sequence may
            be obtained as a Python string with
            "str(parents[i].seq)". Changing this attribute will delete the
            alignment attribute.
        alignment (list[tuple[str, ...]]): Alignment of records, with the
            aligned strings in the same order as the source records. Present if
            and only if the instance is aligned. To set this attribute, call
            the align or set_alignment methods, which will renumber the
            pdb_structure attribute if present.
        p0_aligned (str): records[0] aligned, calculated from alignment.
            Present if and only if instance is aligned.
        pdb_structure (schemarecomb.PDBStructure): Protein Data Bank structure
            that represents the three-dimensional structure of the aligned
            parent sequences.

    Examples:
        To start, you need a FASTA file with at least one parent sequence. For
        these examples, this file is "bgl3_sequence.fasta", which contains the
        HIS-tagged beta-gluocosidase sequence associated with the Protein Data
        Bank (PDB) entry "1GNX".

        The easiest, but slowest way to use this class is to let it handle
        everything through web services. For example, the following script
        will build a six-parent alignment with roughly 70% identity between the
        parents and choose the PDB structure closest to the parents.

        >>> getfixture('bgl3_mock_namespace')
        >>> from schemarecomb import ParentSequences
        >>> fn = 'tests/fixtures/bgl3_1-parent/bgl3_p0.fasta'
        >>> parents = ParentSequences.from_fasta(fn)
        >>> parents.obtain_seqs(6, 0.7)  # BLAST takes about 10 minutes.
        >>> # [sr.name for sr in parents.records]
        >>> parents.align()  # MUSCLE takes about a minute.
        >>> parents.get_PDB()  # BLAST takes about 10 minutes.
        >>> len(parents.records)
        6
        >>> # The following output is shortened for clarity.
        >>> parents.p0_aligned  #doctest: +ELLIPSIS
        '----------------------MHHHHHHMVPAAQQ...WYAEVARTGVLPTA'
        >>> parents.p0_aligned == parents.pdb_structure.renumbering_seq
        True


        You can skip the slow web queries if you already have a FASTA with the
        aligned parents and the PDB structure you want to use:

        >>> from schemarecomb import ParentSequences
        >>> from schemarecomb import PDBStructure
        >>> pdb_fn = 'tests/fixtures/bgl3_full/1GNX.pdb'
        >>> parents_fn = 'tests/fixtures/bgl3_full/bgl3_sequences_aln.fasta'
        >>> pdb = PDBStructure.from_pdb_file(pdb_fn)
        >>> parents = ParentSequences.from_fasta(
        ...     parents_fn,
        ...     pdb_structure=pdb,
        ...     prealigned=True
        ... )
        >>> len(parents.records)
        6
        >>> parents.p0_aligned  #doctest: +ELLIPSIS
        'MHHHHHHMVPAAQQTAMA...RTGVLPTA-----'
        >>> parents.p0_aligned == parents.pdb_structure.renumbering_seq
        True

        You can also save or load a ParentSequences as a JSON:

        >>> from schemarecomb import ParentSequences
        >>> from schemarecomb import PDBStructure
        >>> tempdir = getfixture('tmpdir')  # pytest jargon, ignore this.
        >>> pdb_fn = 'tests/fixtures/bgl3_full/1GNX.pdb'
        >>> parents_fn = 'tests/fixtures/bgl3_full/bgl3_sequences_aln.fasta'
        >>> pdb = PDBStructure.from_pdb_file(pdb_fn)
        >>> parents = ParentSequences.from_fasta(
        ...     parents_fn,
        ...     pdb_structure=pdb,
        ...     prealigned=True
        ... )
        >>> parents_fn = tempdir / 'parents.json'
        >>> parents_json = parents.to_json()
        >>> with open(parents_fn, 'w') as f:
        ...     f.write(parents_json)
        ...
        339501
        >>> with open(parents_fn, 'r') as f:
        ...     parents_json2 = f.read()
        ...
        >>> parents2 = ParentSequences.from_json(parents_json2)
        >>> # parents and parents2 are the same.
        >>> parents.alignment == parents2.alignment
        True

    """

    def __init__(
        self,
        records: list[SeqRecord],
        pdb_structure: Optional['schemarecomb.PDBStructure'] = None,
        auto_align: bool = False,
        prealigned: bool = False
    ) -> None:
        # Input checking.
        if len(records) < 1:
            raise ValueError('records must not be empty.')
        if auto_align and prealigned:
            raise ValueError('auto_align and prealigned must not both be '
                             'True.')

        self.records = records

        self.pdb_structure: 'schemarecomb.PDBStructure'
        if pdb_structure is not None:
            self.pdb_structure = pdb_structure

        if auto_align:
            self.align()
        if prealigned:
            self.new_alignment([str(rec.seq) for rec in records])

    def __setattr__(self, name, value) -> None:
        if name == 'alignment':
            raise AttributeError('The alignment attribute cannot be set '
                                 'directly. Use the align or new_alignment '
                                 'methods.')

        super().__setattr__(name, value)

        if name == 'records' and hasattr(self, '_alignment'):
            # If sequences is changed, alignment is out of date, delete.
            del self._alignment

        if (name == '_alignment' and hasattr(self, 'pdb_structure')) or \
           (name == 'pdb_structure' and hasattr(self, '_alignment')):
            # Renumber existing pdb_structure if alignment changes or if
            # renumber new pdb_structure if self is aligned.
            self.pdb_structure.renumber(self.p0_aligned)

    @property
    def alignment(self) -> list[tuple[str, ...]]:
        try:
            return self._alignment
        except AttributeError:
            raise AttributeError('ParentSequences instance is not aligned.')

    @property
    def p0_aligned(self):
        return ''.join([aminos[0] for aminos in self.alignment])

    def new_alignment(self, aligned_sequences: list[str]) -> None:
        """Add aligned sequences from records to instance.

        Sets the alignment attribute to the input alignments.

        Parameters:
            aligned_sequences: Aligned sequences in the same order as the
                records attribute.

        Raises:
            ValueError: If provided sequences do not match the sequences in the
                records attribute.

        """
        if len(aligned_sequences) != len(self.records):
            raise ValueError('The number of aligned sequences provided must be'
                             ' the same as the number of records.')
        aln_len = len(aligned_sequences[0])
        if any(len(aln_seq) != aln_len for aln_seq in aligned_sequences[1:]):
            raise ValueError('All sequences must be the same length if '
                             'prealigned.')
        for sr, aln_seq in zip(self.records, aligned_sequences):
            if str(sr.seq).replace('-', '') != aln_seq.replace('-', ''):
                raise ValueError('aligned_sequences is does not match '
                                 'self.records.')
        self._alignment = list(zip(*aligned_sequences))

    def align(self, run_locally: bool = False) -> None:
        """Use the MUSCLE web service to align the records.

        Sets the alignment attribute to the resulting alignment. The
        run_locally attribute is currently experimental and probably shouldn't
        be used.

        """
        aln_str = None
        if run_locally:
            try:
                aln_str = local_muscle(self.records)
            except CalledProcessError:
                print('muscle command not available, falling back to MUSCLE '
                      'web service.')

        # Try to use the MUSCLE web API five times.
        for _ in range(5):
            try:
                aln_str = web_muscle(self.records)
                break
            except RuntimeError as e:
                print(e)

        if aln_str is None:
            raise RuntimeError('MUSCLE web API did not work.')

        aln_f = StringIO(aln_str)
        out_SRs = {SR.name: SR for SR in SeqIO.parse(aln_f, 'fasta')}
        reordered_srs = [out_SRs[sr.name] for sr in self.records]

        # This renumbers the self.pdb_structure if present.
        aligned_sequences = [str(sr.seq) for sr in reordered_srs]
        self._alignment = list(zip(*aligned_sequences))

    def obtain_seqs(
        self,
        num_final_sequences: int,
        desired_identity: Optional[float] = None
    ) -> None:
        """Adds new sequences with BLAST.

        Parameters:
            num_final_sequences: Number of desired parent sequences. After
                call, instance should have len(sequences) equal this.
            desired_identity: Desired percent identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.

        Raises:
            ValueError: if num_final_sequences is not greater than
                len(self.records) or if desired_identity is a float and not
                between 0.0 and 1.0, exclusive.

        """
        if len(self.records) >= num_final_sequences:
            len_seqs = len(self.records)
            raise ValueError(f'Requested {num_final_sequences} total sequences'
                             f' but there are already {len_seqs} sequences.')

        # This results in two calls, but we want to validate before BLAST.
        desired_identity = self._check_identity(desired_identity)

        query_sr = self.records[0]
        query_seq = str(query_sr.seq)
        candidate_seqs = list(query_blast(query_seq))
        # query_blast gives iter but list should be small enough

        self.add_from_candidates(candidate_seqs, num_final_sequences,
                                 desired_identity)

    def add_from_candidates(
        self,
        candidate_sequences: list[SeqRecord],
        num_final_sequences: int,
        desired_identity: float = None
    ) -> None:
        """Add new parent sequences from list of candidates.

        Finds the set of sequences in candidate_sequences that gives the
        smallest maximum difference between desired_identity and the calculated
        identity between any two sequences in the set.

        Parameters:
            candidate_sequences: Sequences to choose from.
            num_final_sequences: Number of desired parent sequences. After
                call, instance should have len(sequences) equal this.
            desired_identity: Desired identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.

        Raises:
            ValueError: if num_final_sequences is not greater than
                len(self.records) or if desired_identity is a float and not
                between 0.0 and 1.0, exclusive.

        """

        num_preexisting_parents = len(self.records)
        num_additional = num_final_sequences - num_preexisting_parents
        if num_additional <= 0:
            return
        if len(candidate_sequences) < num_additional:
            raise ValueError('candidate_sequences is not large enough to '
                             f'add {num_additional} sequences.')
        desired_identity = self._check_identity(desired_identity)

        best_cands = choose_candidates(
            candidate_sequences,
            self.records,
            num_additional,
            desired_identity
        )
        self.records += best_cands

    def get_PDB(self) -> None:
        """Construct from ParentSequences using BLAST and PDB.

        The best structure is found by using BLAST to download candidate PDB
        sequences, then the sequence with the largest minimum identity to the
        parents is selected and set to the parent_alignment attribute.

        Parameters:
            parent_aln: Parent alignment used to query the PDB. Note that
                parent_aln[0] is used in the query and PDBStructure alignment.

        Raises:
            ValueError: If no matching PDB structure could be found.

        """
        query_str = str(self.records[0].seq)

        pdb_srs = list(query_blast(query_str, 'pdbaa', 100))

        # Find the PDB struct with largest minimum identity to the parents.
        pdb_min_map = {}
        for pdb_sr in pdb_srs:
            min_iden = min(calc_identity(parent, pdb_sr) for parent
                           in self.records)
            pdb_min_map[pdb_sr.id] = min_iden

        try:
            best_id = max(pdb_min_map, key=lambda pdbid: pdb_min_map[pdbid])
        except ValueError:
            raise RuntimeError('No best PDB could be found.')

        # Get the pdb_structure from rcsb.
        _, acc, chain = best_id.split('|')
        url = 'https://files.rcsb.org/view/' + acc + '.pdb'
        with urlopen(url) as f:
            pdb_structure = schemarecomb.PDBStructure.from_pdb_file(
                f,
                chain=chain
            )

        self.pdb_structure = pdb_structure

    @classmethod
    def from_fasta(cls, fasta_fn: str, **kwargs) -> '_ParentSequences':
        """Contruct instance from FASTA file.

        Parameters:
            fasta_fn: filename of FASTA file, including relative path.
            **kwargs: Additional keyword args for __init__. For example, you
                can specify "auto_align=True" in this constructor.

        Return:
            ParentsSequences instance constructed from input FASTA file.

        """
        seqs = list(SeqIO.parse(fasta_fn, 'fasta'))
        return cls(seqs, **kwargs)

    def to_json(self) -> str:
        """Convert instance to a JSON-formatted string.

        Return:
            Instance converted to a JSON string.

        """

        records = [[str(sr.seq), sr.id, sr.name, sr.description] for sr
                   in self.records]

        # TODO: probably add a version?
        out_list: list = [records]

        try:
            out_list.append(self._alignment)
        except AttributeError:
            out_list.append(None)

        try:
            out_list.append(self.pdb_structure.to_json())
        except AttributeError:
            out_list.append(None)

        return json.dumps(out_list)

    @classmethod
    def from_json(cls, in_json: str) -> '_ParentSequences':
        """Construct instance from JSON.

        Parameters:
            in_json: JSON-formatted string representing a ParentSequences.

        Return:
            ParentSequences instance created from in_json.

        """
        records, alignment, pdb = json.loads(in_json)

        seq_records = []
        for sr_list in records:
            sr_seq, sr_id, sr_name, sr_desc = sr_list
            seq = Seq(sr_seq)
            sr = SeqRecord(seq, id=sr_id, name=sr_name, description=sr_desc)
            seq_records.append(sr)

        if pdb is not None:
            pdb = schemarecomb.PDBStructure.from_json(pdb)

        new_instance = cls(seq_records, pdb)

        if alignment is not None:
            # This will renumber the pdb_structure but that's okay.
            alignment = [tuple(ele) for ele in alignment]
            new_instance._alignment = alignment

        return new_instance

    def _check_identity(self, desired_identity: Optional[float]) -> float:
        """Validate desired_identity, set to records average if None."""
        if desired_identity is None:
            # Follow default behavior: 0.7 if there's only one sequence.
            # Average cross-wise identity of the sequences otherwise.
            if len(self.records) == 1:
                desired_identity = 0.7
            else:
                # calculate average pairwise identity of parents
                identities = [calc_identity(sr1, sr2) for sr1, sr2
                              in combinations(self.records, 2)]
                desired_identity = sum(identities) / len(identities)
        if not 0.0 < desired_identity < 1.0:
            raise ValueError('desired_identity must be between 0.0 and 1.0, '
                             'exclusive.')
        return desired_identity
