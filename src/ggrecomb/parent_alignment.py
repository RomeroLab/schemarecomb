# parent_alignment.py

"""Alignments of parental sequences for combinatorial protein libraries.

This module provides the ParentAlignment class, which represents an alignment
of "parent" protein sequences used in ggrecomb library creation.

Note that the first sequence (ParentAlignment.sequences[0]) has
special importance, as it's used to find additional parental sequences and PDB
structures. For parent sequence sets with high enough identity (~60%), this
should generally not be an issue, but you may get different results by swapping
the order of the sequences.
"""

from collections.abc import Sequence
from itertools import combinations
import json
import os
from subprocess import CalledProcessError, run
from typing import Optional, Union
import urllib.request

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .blast import query_blast
from .pa_candidates import choose_candidates
from ggrecomb import PDBStructure
from .utils import calc_identity


class _ParentAlignment(Sequence):
    """
    Alignment of parent amino acid sequences for recombinant library design.

    This class sets up the data needed to run recombinant design algorithms,
    e.g. intaking parental sequences, finding additional parents, intaking or
    finding PDB strucutures. Instances of this class are passed into functions
    further down the ggrecomb pipeline.

    Parameters:
        parent_SeqRecords: Parental amino acid sequences for ggrecomb
            calculations. Note that the first sequence
            (ParentAlignment.sequences[0]) has special importance, as it's used
            to find additional parental sequences and PDB structures.
        auto_align: If True, the align method is called when the sequences
            attribute changes, including upon initialization.
        prealigned: If True, input sequences are considered aligned upon
            initialization. parent_SeqRecords and aligned_sequences are then
            considered immutable. auto_align and prealigned must not both be
            True.
        pdb_structure: Protein Data Bank structure that represents the three-
            dimensional structure of the aligned parent sequences.

    Attributes:
        sequences: Unaligned parent sequences.
        aligned_sequences: Aligned sequences. Must not be None when instance is
            passed into further methods.
        pdb_structure: Protein Data Bank structure that represents the three-
            dimensional structure of the aligned parent sequences.


    Examples:
        The easiest way to use this class is to initialize via a FASTA file
        with default arguments. This will automatically align the file's
        sequences and the instance will be ready to use downstream:

        >>> from ggrecomb import ParentAlignent
        >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')

        Alternatively, you can use the same procedure and add additional
        sequences with BLAST (might take awhile) by doing:

        >>> from ggrecomb import ParentAlignent
        >>> p_aln.obtain_seqs(num_final_sequences=6, desired_identity=0.65)
        >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')

        You can skip the BLAST search if you have a FASTA file of
        candidate_sequences:

        >>> from Bio import SeqIO
        >>> from ggrecomb import ParentAlignent
        >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')
        >>> candidate_seqs = list(SeqIO.parse('candidates.fasta', 'fasta'))
        >>> p_aln.add_from_candidates(candidate_sequences=candidate_seqs,
        ...                           num_final_sequences=6,
        ...                           desired_identity=0.65)

        Or if you want to build a library from an amino acid sequence stored as
        a string instead:

        >>> from ggrecomb import ParentAlignent
        >>> sequence = 'MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF'
        >>> p_aln = from_single(sequence=sequence,
        ...                     name='YP_025292.1',
        ...                     num_final_sequences=6,
        ...                     desired_identity=0.65)

        You can also save a ParentAlignment as a JSON:

        >>> from ggrecomb import ParentAlignent
        >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')
        >>> aln_json = p_aln.to_json()
        >>> with open('p_aln.json', 'w') as f:
        ...     f.write(aln_json)
        >>> with open('p_aln.json', 'r') as f:
        ...     aln_json2 = f.read()
        >>> p_aln2 = ParentAlignment.from_json(aln_json2)
        >>> # p_aln and p_aln2 are the same
        >>> all(p1_aas == p2_aas for p1_aas, p2_aas in zip(p_aln, p_aln2))
        True

        You can also use the normal sequence operations on the instance
        directly, provided that the instance is aligned. These operations apply
        to the amino acids at each alignment position. For example, this will
        print out the parental amino acids at the  first and tenth positions.

        >>> from ggrecomb import ParentAlignment
        >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')
        >>> print(p_aln[0], p_aln[9])
        ('M', 'M', 'M', 'M', 'M', 'M') ('P', 'A', 'G', 'Y', 'T', 'A')

    """

    def __init__(
        self,
        parent_SeqRecords: list[SeqRecord],
        auto_align: bool = True,
        prealigned: bool = False,
        pdb_structure: 'PDBStructure' = None
    ) -> None:
        # Input checking.
        if len(parent_SeqRecords) < 1:
            raise ValueError('sequences must not be empty.')
        if auto_align and prealigned:
            raise ValueError('auto_align and prealigned must not both be '
                             'True.')
        if prealigned:
            seq_len = len(str(parent_SeqRecords[0].seq))
            if not all(len(str(sr.seq)) == seq_len for sr
                       in parent_SeqRecords[1:]):
                raise ValueError('All sequences must be the same length if '
                                 'prealigned.')

        # Attribute static typing.
        self.aligned_sequences: list[str]
        self.pdb_structure: 'PDBStructure'

        self._auto_align = auto_align
        self._prealigned = prealigned
        self.parent_SeqRecords = parent_SeqRecords
        if pdb_structure is not None:
            self.pdb_structure = pdb_structure

    def __setattr__(self, name, value) -> None:
        super().__setattr__(name, value)

        if name == 'parent_SeqRecords':
            # Prealigned implies immutability.  The exception is when
            # aligned_sequences is not an attr, which should only happen
            # on init.
            if self._prealigned:
                if hasattr(self, 'aligned_sequences'):
                    raise ValueError('sequences cannot be changed if '
                                     'prealigned.')
                else:
                    aln_seqs = [str(sr.seq) for sr in value]
                    self.aligned_sequences = aln_seqs
            else:
                # If sequences is changed, want to modify aligned_sequences.
                if self._auto_align:
                    self.align()
                else:
                    # Delete aligned_sequences and wait for later align() call.
                    del self.aligned_sequences

        elif name == 'aligned_sequences':
            seq_len = len(value[0])
            if not all(len(s) == seq_len for s in value[1:]):
                raise ValueError('All aligned sequences must be the same '
                                 'length.')
            self._amino_acids = tuple(zip(*value))
            if hasattr(self, 'pdb_structure'):
                self.pdb_structure.renumber(self.aligned_sequences[0])

        elif name == 'pdb_structure' and hasattr(self, 'aligned_sequences'):
            # Want to renumber pdb if aligned.
            self.pdb_structure.renumber(self.aligned_sequences[0])

    def __delattr__(self, name):
        super().__delattr__(name)
        if name == 'aligned_sequences':
            del self._amino_acids
            if hasattr(self, 'pdb_structure'):
                self.pdb_structure.derenumber()

    def __getitem__(self, key):
        """Amino acids at position key in the alignment.

        Raises:
            AttributeError if the instance has not been aligned.
        """
        try:
            return self._amino_acids[key]
        except AttributeError:
            raise AttributeError('ParentAlignment is not aligned yet.')

    def __len__(self):
        """Number of amino acid positions in the alignment.

        Raises AttributeError if the instance has not been aligned.
        """
        try:
            return len(self._amino_acids)
        except AttributeError:
            raise AttributeError('ParentAlignment is not aligned yet.')

    @classmethod
    def from_fasta(cls, fasta_fn: str, **kwargs) -> '_ParentAlignment':
        """Contruct instance from FASTA file.

        Parameters:
            fasta_fn: filename of FASTA file, including relative path.
            **kwargs: Additional keyword args for __init__. For example, you
                can specify "auto_align=False" in this constructor.

        """
        seqs = list(SeqIO.parse(fasta_fn, 'fasta'))
        return cls(seqs, **kwargs)

    @classmethod
    def from_single(cls,
                    sequence: Union[str, SeqRecord],
                    name: Optional[str] = None,
                    num_final_sequences: int = 3,
                    desired_identity: float = 0.7,
                    **kwargs) -> '_ParentAlignment':
        """Contruct instance from single amino acid sequence.

        Parameters:
            sequence: First parent sequence. Will be used to find additional
                sequences and PDB structure.
            name: Name of sequence. Ignored if sequence is a SeqRecord; must be
                specified if sequence is a string.
            num_final_sequences: Number of desired parent sequences. Returned
                ParentAlignment instance should have len(sequences) equal this.
            **kwargs: Additional keyword args for __init__. For example, you
                can specify "auto_align=False" in this constructor.

        """
        if isinstance(sequence, str):
            if name is None:
                raise ValueError('name parameter is required when sequence is '
                                 'a string.')
            sequence = SeqRecord(Seq(sequence), id=name, name=name)
        pa = cls([sequence], **kwargs)
        pa.obtain_seqs(num_final_sequences, desired_identity)
        return pa

    def obtain_seqs(self,
                    num_final_sequences: int,
                    desired_identity: Optional[float] = None) -> None:
        """Adds new sequences with BLAST.

        Parameters:
            num_final_sequences: Number of desired parent sequences. After
                call, instance should have len(sequences) equal this.
            desired_identity: Desired percent identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.
        """
        if len(self.parent_SeqRecords) >= num_final_sequences:
            len_seqs = len(self.parent_SeqRecords)
            raise ValueError(f'Requested {num_final_sequences} total sequences'
                             f' but there are already {len_seqs} sequences.')

        # This results in two calls, but we want to validate before BLAST.
        desired_identity = self._check_identity(desired_identity)

        query_sr = self.parent_SeqRecords[0]
        query_seq = str(query_sr.seq)
        candidate_seqs = list(query_blast(query_seq))
        # blast_query gives iter but list should be small enough

        self.add_from_candidates(candidate_seqs, num_final_sequences,
                                 desired_identity)

    def add_from_candidates(self,
                            candidate_sequences: list[SeqRecord],
                            num_final_sequences: int,
                            desired_identity: float = None) -> None:
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
        """

        num_preexisting_parents = len(self.parent_SeqRecords)
        num_additional = num_final_sequences - num_preexisting_parents
        if num_additional <= 0:
            return
        if len(candidate_sequences) < num_additional:
            raise ValueError('candidate_sequences is not larget enough to '
                             f'add {num_additional} sequences.')
        desired_identity = self._check_identity(desired_identity)

        best_cands = choose_candidates(candidate_sequences,
                                       self.parent_SeqRecords,
                                       num_additional, desired_identity)
        self.parent_SeqRecords += best_cands

    def align(self) -> None:
        """Align sequences with MUSCLE and set to aligned_sequences."""
        # TODO: add MUSCLE web server fallback.
        print('running alignment')

        IN_FN = 'temp_muscle_input.fasta'
        OUT_FN = 'temp_muscle_output.fasta'

        # make file for MUSCLE input
        SeqIO.write(self.parent_SeqRecords, IN_FN, 'fasta')
        in_labels = [SR.name for SR in SeqIO.parse(IN_FN, 'fasta')]

        # run MUSCLE
        try:
            run(f'muscle -in {IN_FN} -out {OUT_FN}', shell=True, check=True)
        except CalledProcessError:
            print('Something is wrong with MUSCLE call. Is MUSCLE installed?')
            raise
        finally:
            os.remove(IN_FN)

        # need to reorder MUSCLE output to be consistent with input
        out_SRs = {SR.name: SR for SR in SeqIO.parse(OUT_FN, 'fasta')}
        reordered_srs = [out_SRs[label] for label in in_labels]
        self.aligned_sequences = [str(sr.seq) for sr in reordered_srs]

        print('Muscle done, cleaning up files')
        os.remove(OUT_FN)

        # Renumber pdb_structure if available.
        if hasattr(self, 'pdb_structure') and self.pdb_structure is not None:
            self.pdb_structure.renumber(self.aligned_sequences[0])

    def get_PDB(self):
        """Construct from ParentAlignment using BLAST and PDB.

        The best structure is found by using BLAST to download candidate PDB
        sequences, then the sequence with the largest minimum identity to the
        parents is selected. This sequence is aligned with the first parent to
        renumber the PDB residues to parent alignment.

        Note that the current implementation doesn't renumber PDB residues that
        aren't found in the first parent. This effectively means that the first
        parent should always be selected from the PDB database.

        Parameters:
            parent_aln: Parent alignment used to query the PDB. Note that
                parent_aln[0] is used in the query and PDBStructure alignment.

        """
        try:
            query_str = self.aligned_sequences[0]
        except AttributeError:
            raise AttributeError('ParentAlignment is not aligned.')

        pdb_srs = list(query_blast(query_str, 'pdbaa', 100))

        # Find the PDB struct with largest minimum identity to the parents.
        best_id = None  # track the best sequence
        best_min_iden = 0.0  # minimum parental identity of best sequence
        for pdb_sr in pdb_srs:
            min_iden = 1.0  # minimum parental identity of current sequence
            is_minimum = True  # whether the current sequence is the best
            for par_sr in self.parent_SeqRecords:
                iden = calc_identity(pdb_sr, par_sr)
                if iden < best_min_iden:  # sequence suboptimal, short-circuit
                    is_minimum = False
                    break
                min_iden = min(min_iden, iden)
            if is_minimum:  # never set false, must be optimal so far
                best_id = pdb_sr.id
                best_min_iden = min_iden

        if best_id is None:
            raise ValueError('No best PDB found.')

        # Get the pdb_structure from rcsb.
        acc = best_id.split('|')[1]
        url = 'https://files.rcsb.org/view/' + acc + '.pdb'
        request = urllib.request.Request(url)
        with urllib.request.urlopen(request) as f:
            pdb_structure = PDBStructure.from_pdb_file(f)

        # pdb_structure.renumber(query_str)  # probably don't want this?

        return pdb_structure

    def to_json(self) -> str:
        """Convert instance to JSON."""
        '''
        if self.aligned_sequences is None:
            seq_records = self.parent_SeqRecords
            aligned = False
        else:
            seq_records = self.aligned_sequences
            aligned = True
        '''

        # TODO: probably add a version?
        unaln_sr_jsons = []
        for sr in self.parent_SeqRecords:
            sr_dict = {'seq': str(sr.seq), 'id': sr.id, 'name': sr.name,
                       'description': sr.description}
            unaln_sr_jsons.append(sr_dict)

        out_dict = {
            'parent_SeqRecords': unaln_sr_jsons,
            'auto_align': self._auto_align,
            'prealigned': self._prealigned
        }

        if hasattr(self, 'aligned_sequences'):
            out_dict['aligned_sequences'] = self.aligned_sequences

        if hasattr(self, 'pdb_structure'):
            out_dict['pdb_structure'] = self.pdb_structure.to_json()

        return json.dumps(out_dict)

    @classmethod
    def from_json(cls, in_json: str) -> '_ParentAlignment':
        """Construct instance from JSON."""
        # TODO: Rewrite these.
        in_dict = json.loads(in_json)

        seq_records = []
        for sr_dict in in_dict['seq_records']:
            seq = Seq(sr_dict['seq'])
            sr_id = sr_dict['id']
            sr_name = sr_dict['name']
            sr_desc = sr_dict['description']
            sr = SeqRecord(seq, id=sr_id, name=sr_name, description=sr_desc)
            seq_records.append(sr)

        # new ParentAlignment with auto_align temporarily disabled
        new_instance = cls(seq_records, auto_align=False)

        if in_dict['aligned']:
            # Move sequences to aligned_sequences and recalculate unaligned.
            # This is not the same thing as prealigned.
            srs = []
            for sr in seq_records:
                seq = Seq(str(sr.seq).replace('-', ''))
                sr = SeqRecord(seq, id=sr.id, name=sr.name,
                               description=sr.description)
                srs.append(sr)
            new_instance.parent_SeqRecords = srs
            new_instance.aligned_sequences = [str(sr.seq) for sr in srs]

        new_instance._auto_align = in_dict['auto_align']
        return new_instance

    def _check_identity(self, desired_identity: Optional[float]) -> float:
        """Validate the input desired_identity, set to default if None."""
        if desired_identity is None:
            # Follow default behavior: 0.7 if there's only one sequence.
            # Average cross-wise identity of the sequences otherwise.
            if len(self.parent_SeqRecords) == 1:
                desired_identity = 0.7
            else:
                # calculate average pairwise identity of parents
                identities = [calc_identity(sr1, sr2) for sr1, sr2
                              in combinations(self.parent_SeqRecords, 2)]
                desired_identity = sum(identities) / len(identities)
        if not 0.0 < desired_identity < 1.0:
            raise ValueError('desired_identity must be between 0.0 and 1.0, '
                             'exclusive.')
        return desired_identity
