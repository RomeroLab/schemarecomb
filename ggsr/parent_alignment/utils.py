# utils.py
"""Utilities module for parent_alignment."""

from Bio import pairwise2, SeqRecord


def _calc_identity(sr1: SeqRecord.SeqRecord, sr2: SeqRecord.SeqRecord) \
        -> float:
    """Calculate the BLAST identity between two sequences.

    Note that "BLAST identity" is defined as the number of matches in a
    pairwise alignment divided by the length of alignment. This is different
    from the "traditional" definition of identity: number of matches divided by
    the minimum of the sequence lengths.
    """
    s1 = str(sr1.seq)
    s2 = str(sr2.seq)
    aln_score = pairwise2.align.globalxx(s1, s2, score_only=True)
    aln_len = len(s1) + len(s2) - aln_score
    return aln_score / aln_len


def iden_diff(sr1, sr2, target_identity):
    return abs(_calc_identity(sr1, sr2) - target_identity)
