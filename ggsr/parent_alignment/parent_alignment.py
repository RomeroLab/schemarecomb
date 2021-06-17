# parent_alignment.py

"""Alignments of parental sequences for combinatorial protein libraries."""
# TODO: MUSCLE

from Bio import Seq, SeqIO, SeqRecord
from Bio import pairwise2


class ParentAlignment:
    def __init__(self, sequences, auto_align=False):
        assert len(sequences) > 0
        self.auto_align = auto_align
        self.sequences = sequences
        self.aligned_sequences = None  # whether muscle has been run

    def __setattr__(self, name, value):
        super().__setattr__(name, value)

        # If sequences is changed, want to modify aligned_sequences.
        if name == 'sequences':
            if self.auto_align:
                self.align()
            else:
                self.aligned_sequences = None

    @classmethod
    def from_fasta(cls, fasta_fn, **args):
        seqs = list(SeqIO.parse(fasta_fn, 'fasta'))
        return cls(seqs, **args)

    @classmethod
    def from_single(cls, sequence, name=None, file=False, num_sequences=3,
                    desired_identity=0.7, **args):
        """Makes alignment from single sequence."""
        if isinstance(sequence, str):
            if name is None:
                raise ValueError('name parameter is required when sequence is '
                                 'a string.')
            sequence = SeqRecord.SeqRecord(Seq.Seq(sequence), id=name)
        elif isinstance(sequence, SeqRecord):
            if name is not None:
                sequence.id = name
        else:
            raise TypeError('sequence parameter must be a String or '
                            'Bio.SeqRecord.')
        pa = cls(sequence, **args)
        pa.obtain_seqs(num_sequences, desired_identity)
        return pa

    def obtain_seqs(self, num_final_sequences, desired_identity=None):
        """Adds new sequences with BLAST.

        Args:
            num_final_sequences: Number of desired parent sequences.
            desired_identity: Desired percent identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.
        """
        from .blast_query import blast_query
        num_preexisting_parents = len(self.sequences)
        num_additional = num_final_sequences - num_preexisting_parents
        if num_additional <= 0:
            return
        if desired_identity is None:
            if num_preexisting_parents == 1:
                desired_identity = 0.7
            else:
                desired_identity = _average_identity(self.sequences)
        else:
            assert 0.0 < desired_identity < 1.0

        query_seq = self.sequences[0]
        candidate_seqs_iter = blast_query(query_seq)

        self.add_from_candidates(candidate_seqs_iter, num_final_sequences,
                                 desired_identity)

    def add_from_candidates(self, candidate_sequences, num_final_sequences,
                            desired_identity=None):
        """Add new parent sequences from list of candidates.

        Finds the set of sequences in candidate_sequences that gives the
        smallest maximum difference between desired_identity and the calculated
        identity between any two sequences in the set.

        Args:
            candidate_sequences: Iterable collection of sequences to choose
                from.
            num_final_sequences: Number of desired parent sequences.
            desired_identity: Desired identity of new sequences
                compared to first parent sequence. Default: if there are
                multiple parent sequences already, this value is taken as the
                average identity between the first parent and other sequences.
                Otherwise, 70% identity.
        """
        # TODO: candidate_sequences should be iterable of SeqRecords, not
        # str (probably a concern with whole class)

        # TODO: repeat code, can probably refactor
        num_preexisting_parents = len(self.sequences)
        num_additional = num_final_sequences - num_preexisting_parents
        if num_additional <= 0:
            return
        if desired_identity is None:
            if num_preexisting_parents == 1:
                desired_identity = 0.7
            else:
                desired_identity = _average_identity(self.sequences)

        # from collections.abc import Iterator
        # if not isinstance(candidate_sequences, Iterator):
        #     candidate_sequences = iter(candidate_sequences)

        import pickle
        piden = desired_identity
        cand_diffs = []
        print('calcing')
        for i, cand in enumerate(candidate_sequences):
            print(i, '\r', end='')
            max_diff = max(abs(piden - _calc_identity(cand, x))
                           for x in self.sequences)
            cand_diffs.append((cand, max_diff))

        # sorted is fast enough (but could be faster with heapq)
        sorted_cand_diffs = list(sorted(cand_diffs, key=lambda x: x[1]))
        print('dumping')
        pickle.dump(sorted_cand_diffs,
                    open('../../tests/bgl3_sample/sorted_cand_diffs.pkl',
                         'wb'))

        print('loading')
        sorted_cand_diffs = pickle.load(
                    open('../../tests/bgl3_sample/sorted_cand_diffs.pkl',
                         'rb'))
        print('done loading')

        from choose_candidates import choose_candidates
        best_cands = choose_candidates(sorted_cand_diffs, num_additional)

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

        print(len(self.sequences))
        self.sequences += best_cands
        print(len(self.sequences))

        self.aligned = False  # new sequences added, need to run muscle again

    def align(self):
        print('running alignment')
        import os
        from subprocess import CalledProcessError, run

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


def _calc_identity(sr1, sr2):
    s1 = str(sr1.seq)
    s2 = str(sr2.seq)
    aln_result2 = pairwise2.align.globalxx(s1, s2, score_only=True)
    return aln_result2 / (len(s1) + len(s2) - aln_result2)


def _average_identity(sequences):
    """Average identity between sequences.

    Helper method for ParentAlignment.obtain_seqs.
    """
    identities = []
    for i, seq1 in enumerate(sequences):
        for seq2 in sequences[i+1:]:
            identities.append(_calc_identity(seq1, seq2))
    return sum(identities) / len(identities)


if __name__ == '__main__':
    aln = ParentAlignment.from_fasta('../../tests/bgl3_sample/bgl3_p1-2.fasta',
                                     auto_align=True)
    cands = list(SeqIO.parse('../../tests/bgl3_sample/bgl3_p3-6.fasta',
                             'fasta'))
    aln.add_from_candidates(cands, 4)
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
