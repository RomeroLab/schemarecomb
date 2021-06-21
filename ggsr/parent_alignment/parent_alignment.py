# parent_alignment.py

"""Alignments of parental sequences for combinatorial protein libraries.

This module provides the ParentAlignment class, which represents an alignment
of "parent" protein sequences used in Golden Gate SCHEMA-RASPPi (GGSR) library
creation.

Note that the first sequence (ParentAlignment.sequences[0]) has
special importance, as it's used to find additional parental sequences and PDB
structures. For parent sequence sets with high enough identity (~60%), this
should generally not be an issue, but you may get different results by swapping
the order of the sequences.
"""

from itertools import combinations
import json
import os
from subprocess import CalledProcessError, run
from typing import Optional, Union

from Bio import Seq, SeqIO, SeqRecord

from .blast_query import blast_query
from .choose_candidates import choose_candidates
from .utils import _calc_identity


class ParentAlignment:
    """Alignment of parental protein sequences for GGSR.

    Public attributes:
        sequences: Unaligned parent sequences.
        aligned_sequences: Aligned sequences. Will be None if
            auto_align==False, or automatically assigned if auto_align=True,
            each time sequences is set. Must not be None when instance is
            passed into further methods.
        auto_align: Whether align() is called or aligned_sequences is set to
            None when setting sequences. Increase performance if set to False.

    Public methods:
        obtain_seqs: BLAST search to find and add candidate sequences.
        add_from_candidates: Choose to add from a list of candidate sequences.
        align: Align sequences with MUSCLE and set to aligned_sequences.
        to_json: Convert instance to JSON.

    Constructors:
        from_fasta: Construct instance from FASTA file.
        from_single: Construct instance from sequence string or SeqRecord.
        from_json: Construct instance from JSON.

    Instances of this class are passed into functions further down the GGSR
    pipeline.

    Note that the first sequence (ParentAlignment.sequences[0]) has special
    importance, as it's used to find additional parental sequences and PDB
    structures.

    The easiest way to use this class is to initialize via a FASTA file with
    default arguments. This will automatically align the file's sequences and
    the instance will be ready to use downstream:
    >>> from ggsr.parent_alignment import ParentAlignent
    >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')

    Alternatively, you can use the same procedure and add additional sequences
    with BLAST (might take awhile) by doing:
    >>> from ggsr.parent_alignment import ParentAlignent
    >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')
    >>> p_aln.obtain_seqs(num_final_sequences=6, desired_identity=0.65)

    You can skip the BLAST search if you have a FASTA file of
    candidate_sequences:
    >>> from Bio import SeqIO
    >>> from ggsr.parent_alignment import ParentAlignent
    >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')
    >>> candidate_seqs = list(SeqIO.parse('candidates.fasta', 'fasta'))
    >>> p_aln.add_from_candidates(candidate_sequences=candidate_seqs,
    ...                           num_final_sequences=6,
    ...                           desired_identity: 0.65)

    Or if you want to build a library from an amino acid sequence stored as a
    string instead:
    >>> from ggsr.parent_alignment import ParentAlignent
    >>> sequence = 'MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF'
    >>> p_aln = from_single(cls,
    ...                     sequence=sequence,
    ...                     name='YP_025292.1',
    ...                     num_final_sequences=6,
    ...                     desired_identity=0.65)

    You can also save a ParentAlignment as a JSON:
    >>> from ggsr.parent_alignment import ParentAlignent
    >>> p_aln = ParentAlignment.from_fasta('P450_AA_sequences.fasta')
    >>> aln_json = p_aln.to_json()
    >>> with open('p_aln.json', 'w') as f:
    ...     f.write(aln_json)
    >>> with open('p_aln.json', 'r') as f:
    ...     aln_json2 = f.read()
    >>> p_aln2 = ParentAlignment.from_json(aln_json2)
    >>> # p_aln and p_aln2 should be the same
    """

    def __init__(self, sequences: list[SeqRecord.SeqRecord],
                 auto_align: bool = True) -> None:
        """Initalize ParentAlignment.

        Args:
            sequences: Parental amino acid sequences for GGSR calculations.
            auto_align: Whether align is called when sequences is reset. If
                False, aligned_sequences is set to None when sequences is set.
                Note that align will be called upon initialization if
                auto_align == True.
        """
        assert len(sequences) > 0
        self.auto_align = auto_align
        self.aligned_sequences = None  # not needed here, included for clarity
        self.sequences = sequences  # must be set last since it affects others

    def __setattr__(self, name, value) -> None:
        """Set aligned_sequences according to auto_align if sequences reset."""
        super().__setattr__(name, value)

        # If sequences is changed, want to modify aligned_sequences.
        if name == 'sequences':
            if self.auto_align:
                self.align()
            else:
                self.aligned_sequences = None

    @classmethod
    def from_fasta(cls, fasta_fn: str, **kwargs) -> 'ParentAlignment':
        """Contruct instance from FASTA file.

        Args:
            fasta_fn: filename of FASTA file, including relative path.
            **kwargs: Additional keyword args for __init__. For example, you
                can specify "auto_align=False" in this constructor.
        """
        seqs = list(SeqIO.parse(fasta_fn, 'fasta'))
        return cls(seqs, **kwargs)

    @classmethod
    def from_single(cls,
                    sequence: Union[str, SeqRecord.SeqRecord],
                    name: Optional[str] = None,
                    num_final_sequences: int = 3,
                    desired_identity: float = 0.7,
                    **kwargs) -> 'ParentAlignment':
        """Contruct instance from single amino acid sequence.

        Args:
            sequence: First parent sequence. Will be used to find additional
                sequences and PDB structure.
            name: Name of sequence. May be None if sequence is a SeqRecord;
                must be specified if sequence is a string.
            num_final_sequences: Number of desired parent sequences. Returned
                ParentAlignment instance should have len(sequences) equal this.
            **kwargs: Additional keyword args for __init__. For example, you
                can specify "auto_align=False" in this constructor.
        """
        if isinstance(sequence, str):
            if name is None:
                raise ValueError('name parameter is required when sequence is '
                                 'a string.')
            sequence = SeqRecord.SeqRecord(Seq.Seq(sequence), id=name)
        elif isinstance(sequence, SeqRecord.SeqRecord):
            if name is not None:
                sequence.id = name
        else:
            raise TypeError('sequence parameter must be a String or '
                            'Bio.SeqRecord.SeqRecord.')
        pa = cls(sequence, **kwargs)
        pa.obtain_seqs(num_final_sequences, desired_identity)
        return pa

    def obtain_seqs(self,
                    num_final_sequences,
                    desired_identity: Optional[float] = None) -> None:
        """Adds new sequences with BLAST.

        Args:
            num_final_sequences: Number of desired parent sequences. After
                call, instance should have len(sequences) equal this.
            desired_identity: Desired percent identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.
        """
        query_seq = self.sequences[0]
        candidate_seqs_iter = blast_query(query_seq)

        self.add_from_candidates(candidate_seqs_iter, num_final_sequences,
                                 desired_identity)

    def add_from_candidates(self,
                            candidate_sequences: list[SeqRecord],
                            num_final_sequences: int,
                            desired_identity: float = None) -> None:
        """Add new parent sequences from list of candidates.

        Finds the set of sequences in candidate_sequences that gives the
        smallest maximum difference between desired_identity and the calculated
        identity between any two sequences in the set.

        Args:
            candidate_sequences: Sequences to choose from.
            num_final_sequences: Number of desired parent sequences. After
                call, instance should have len(sequences) equal this.
            desired_identity: Desired identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.
        """

        num_preexisting_parents = len(self.sequences)
        num_additional = num_final_sequences - num_preexisting_parents
        if num_additional <= 0:
            return
        if desired_identity is None:
            if num_preexisting_parents == 1:
                desired_identity = 0.7
            else:
                # calculate average pairwise identity of parents
                identities = [_calc_identity(sr1, sr2) for sr1, sr2
                              in combinations(self.sequences, 2)]
                desired_identity = sum(identities) / len(identities)
        assert 0.0 < desired_identity < 1.0

        best_cands = choose_candidates(candidate_sequences, self.sequences,
                                       num_additional, desired_identity)

        # TODO: remove this
        # naive way, probably a good test?
        '''
        from itertools import combinations
        cands = sorted_cand_diffs[:20]
        best_diff = 1.0
        best_set = None
        print(1140)
        for k, srs in enumerate(combinations(cands, 3)):
            print(k, '\r', end='')
            # max_diff = 0.0
            max_diff = max([sr[1] for sr in srs])
            for i, (sr, sra) in enumerate(srs):
                # for p in self.sequences:
                #     diff = abs(_calc_identity(sr, p) - 0.7)
                #     if diff > max_diff:
                #         max_diff = diff
                for sr2, _ in srs[i+1:]:
                    diff = abs(_calc_identity(sr, sr2) - 0.7)
                    if diff > max_diff:
                        max_diff = diff
            if max_diff < best_diff:
                best_diff = max_diff
                best_set = srs
        print([bs[0].id for bs in best_set])
        '''

        self.sequences += best_cands

    def align(self) -> None:
        """Align sequences with MUSCLE and set to aligned_sequences."""
        print('running alignment')

        IN_FN = 'temp_muscle_input.fasta'
        OUT_FN = 'temp_muscle_output.fasta'

        # make file for MUSCLE input
        SeqIO.write(self.sequences, IN_FN, 'fasta')
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
        self.aligned_sequences = [out_SRs[label] for label in in_labels]

        print('Muscle done, cleaning up files')
        os.remove(OUT_FN)

    def to_json(self) -> str:
        """Convert instance to JSON."""
        if self.aligned_sequences is None:
            seq_records = self.sequences
            aligned = False
        else:
            seq_records = self.aligned_sequences
            aligned = True

        # TODO: probably add a version?
        out_dict = {'seq_records': [], 'aligned': aligned,
                    'auto_align': self.auto_align}
        # important to preserve order of sequences
        for sr in seq_records:
            sr_dict = {'seq': str(sr.seq), 'id': sr.id, 'name': sr.name,
                       'description': sr.description}
            out_dict['seq_records'].append(sr_dict)
        return json.dumps(out_dict)

    @classmethod
    def from_json(cls, in_json: str) -> 'ParentAlignment':
        """Construct instance from JSON."""
        in_dict = json.loads(in_json)

        seq_records = []
        for sr_dict in in_dict['seq_records']:
            seq = Seq.Seq(sr_dict['seq'])
            sr_id = sr_dict['id']
            sr_name = sr_dict['name']
            sr_desc = sr_dict['description']
            sr = SeqRecord.SeqRecord(seq, id=sr_id, name=sr_name,
                                     description=sr_desc)
            seq_records.append(sr)

        # new ParentAlignment with auto_align temporarily disabled
        new_instance = cls(seq_records, auto_align=False)

        if in_dict['aligned']:
            # move sequences to aligned_sequences and recalculate unaligned
            unaligned_seqs = []
            for sr in seq_records:
                u_str = str(sr.seq).replace('-', '')
                u_seq = Seq.Seq(u_str)
                u_sr = SeqRecord.SeqRecord(u_seq, id=sr.id, name=sr.name,
                                           description=sr.description)
                unaligned_seqs.append(u_sr)
            new_instance.sequences = unaligned_seqs
            new_instance.aligned_sequences = seq_records

        new_instance.auto_align = in_dict['auto_align']
        return new_instance


if __name__ == '__main__':
    aln = ParentAlignment.from_fasta('../../tests/bgl3_sample/bgl3_p1-2.fasta',
                                     auto_align=True)
    aln_json = aln.to_json()
    aln2 = ParentAlignment.from_json(aln_json)
    print(aln)
    print(aln2)

    # cands = list(SeqIO.parse('../../tests/bgl3_sample/bgl3_p3-6.fasta',
    #                          'fasta'))
    # aln.add_from_candidates(cands, 4)
    # aln.align()
    '''
    import blast_query
    qs = list(SeqIO.parse('../../tests/bgl3_sample/bgl3_p1-2.fasta',
                          'fasta'))[0].seq
    candidate_srs = blast_query.blast_query(qs)
    print('writing')
    SeqIO.write(candidate_srs, '../../tests/bgl3_sample/query_seqs2.fasta',
                'fasta')
    '''
    '''
    cands = list(SeqIO.parse('../../tests/bgl3_sample/query_seqs.fasta',
                             'fasta'))
    aln.add_from_candidates(cands, 5, 0.7)
    '''
