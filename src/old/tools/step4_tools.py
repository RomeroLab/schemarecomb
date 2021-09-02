from collections import namedtuple
import itertools
import re

from Bio import Seq, SeqRecord

from tools import general_tools


def first_overhang_seq(overhang, str_seq):
    """ confusing, but the first overhang will overlap with the first codon to
    some degree. Since the first breakpoint is at 0, AA2 for this breakpoint is
    at sequence position 0, while AA1 is at an imaginary sequence position -1.
    oh_pos represents the overhang position relative to AA1, so the number of
    overlapping positions is oh_in_seq=oh_pos+1. Therefore, the last oh_in_seq
    bases in the overhang and the first oh_in_seq bases in the sequence should
    be the same. It follows that we will need to add the 4-oh_in_seq leftmost
    overhang bases to the final fragment.
    """
    oh_pos, oh_seq = overhang
    if oh_pos == -1:
        return oh_seq  # just append oh onto sequence
    # number of bases on left side of seq that should be on right side of oh
    oh_in_seq = oh_pos + 1
    assert oh_seq[-oh_in_seq:] == str_seq[:oh_in_seq]
    return oh_seq[:4-oh_in_seq]


def last_overhang_seq(overhang, str_seq):
    """ confusing, but the last overhang will overlap with the last codon to
    some degree. Since the last breakpoint is at len(alignment), AA1 for this
    breakpoint is at the last sequence position, while AA2 is at an imaginary
    sequence position len(alignment). oh_pos represents the overhang position
    relative to AA1, so the number of overlapping positions is
    oh_in_seq=3-oh_pos. Therefore, the first oh_in_seq bases in the overhang
    and the last oh_in_seq bases in the sequence should be the same. It follows
    that we will need to add the 4-oh_in_seq rightmost overhang bases to the
    final fragment.
    """
    oh_pos, oh_seq = overhang
    if oh_pos == 3:
        return oh_seq  # just append oh onto sequence
    # number of bases on right side of seq that should be in left side of oh
    oh_in_seq = 3 - oh_pos
    assert oh_seq[:oh_in_seq] == str_seq[-oh_in_seq:]
    return oh_seq[oh_in_seq-4:]


def _get_valid_codon(AA, pattern):
    residue_cdns = general_tools.rev_code[AA]
    possible_cdns = [c for c in residue_cdns if re.fullmatch(pattern, c)]
    cdn = possible_cdns[0]
    return cdn


def overhang_CDN_seq(overhang, cdn1, cdn2, block_front_back):
    position, oh_seq = overhang

    # build codon patterns for both codons
    overall_pattern = '.' * position + oh_seq + '.' * (2 - position)
    pattern1 = overall_pattern[:3]
    pattern2 = overall_pattern[3:]

    AA1 = general_tools.code[cdn1]
    AA2 = general_tools.code[cdn2]

    cdn1 = _get_valid_codon(AA1, pattern1)
    cdn2 = _get_valid_codon(AA2, pattern2)

    cdn_seq = cdn1 + cdn2

    assert block_front_back in ('front', 'back')
    if block_front_back == 'front':
        return cdn_seq[position:]

    if position == 2:
        return cdn_seq
    return cdn_seq[:-2 + position]


# Script for confirming order against amino acid sequence below

class Restriction_Enzyme(object):
    def __init__(self, rec_site, offset, cut_len, sticky_end_strand):
        self.site = rec_site
        self.offset = offset
        self.cut_len = cut_len

        # which strand sticky end stays on relative to site, 'same' or 'rc'
        assert sticky_end_strand == 'same' or sticky_end_strand == 'rc'
        self.sticky_end_strand = sticky_end_strand


bsaI = Restriction_Enzyme('GGTCTC', 7, 4, 'rc')

Overhang_set = namedtuple('Overhang_set', ['fwd_start', 'fwd_end', 'rev_start',
                                           'rev_end'])


def same_strand_cut(seq, rs_site, enzyme):
    se_len = enzyme.cut_len  # length of sticky ends
    offset = enzyme.offset
    if enzyme.sticky_end_strand == 'same':
        cut_site = rs_site + offset + se_len
        new_seq = seq[cut_site:]
    else:
        cut_site = rs_site + offset
        cut_seq = seq[cut_site:]
        new_seq = cut_seq[:se_len].lower() + cut_seq[se_len:]
    return new_seq


def rc_strand_cut(seq, rs_site, enzyme):
    # assume same and rc srands end at same position
    se_len = enzyme.cut_len  # length of sticky ends
    offset = enzyme.offset
    if enzyme.sticky_end_strand == 'same':
        cut_site = -rs_site - offset
        cut_seq = seq[:cut_site]
        new_seq = cut_seq[:-se_len] + cut_seq[-se_len:].lower()
    else:
        cut_site = -rs_site - offset - se_len
        new_seq = seq[:cut_site]
    return new_seq


def reverse_complements(oh1, oh2):
    s1 = Seq.Seq(oh1)
    s2 = Seq.Seq(oh2)
    return s1 == s2.reverse_complement()


def concatenate_frags(frag1, frag2):
    # order matters, result is frag1 + frag2
    frag1_overhangs = frag1.get_overhangs()
    frag2_overhangs = frag2.get_overhangs()
    assert frag1_overhangs.fwd_end or frag2_overhangs.rev_start
    if frag1_overhangs.fwd_end:
        frag1_oh = frag1_overhangs.fwd_end
        frag2_oh = frag2_overhangs.rev_end
        assert reverse_complements(frag1_oh, frag2_oh)
        oh_len = len(frag1_oh)
        fwd_seq = frag1.fwd[:-oh_len] + frag1_oh.upper() + frag2.fwd
        rev_seq = frag2.rev[:-oh_len] + frag2_oh.upper() + frag1.rev
    else:
        frag1_oh = frag1_overhangs.rev_start
        frag2_oh = frag2_overhangs.fwd_start
        assert reverse_complements(frag1_oh, frag2_oh)
        oh_len = len(frag1_oh)
        fwd_seq = frag1.fwd + frag2_oh.upper() + frag2.fwd[oh_len:]
        rev_seq = frag2.rev + frag1_oh.upper() + frag1.rev[oh_len:]

    return Fragment(fwd_seq, rev_seq)


def get_seq_overhangs(seq):
    start_se_len = re.search('[A-Z]', seq).start()  # fwd seq beginning
    end_se_len = re.search('[A-Z]', seq[::-1]).start()  # fwd seq end
    if end_se_len:
        return seq[:start_se_len], seq[-end_se_len:]
    return seq[:start_se_len], ''


class Fragment(object):
    def __init__(self, fwd_seq: str, rev_seq: str):
        self.fwd = fwd_seq
        self.rev = rev_seq

    @classmethod
    def from_SeqRecord(cls, sr: SeqRecord.SeqRecord) -> 'Fragment':
        fwd = str(sr.seq).upper()
        rev = str(sr.seq.reverse_complement()).upper()
        return cls(fwd, rev)

    def digest(self, enzyme: Restriction_Enzyme):
        # find restriction sites
        fwd_sites = [m.start() for m in re.finditer('(?=%s)' % enzyme.site,
                                                    self.fwd)]
        if len(fwd_sites) > 1:
            raise Exception('More than one foward restriction site found.')
        fwd_site_index = fwd_sites[0]
        rev_sites = [m.start() for m in re.finditer('(?=%s)' % enzyme.site,
                                                    self.rev)]
        if len(rev_sites) > 1:
            raise Exception('More than one reverse restriction site found.')
        rev_site_index = rev_sites[0]

        # create new forward seq
        new_fwd = self.fwd
        new_fwd = same_strand_cut(new_fwd, fwd_site_index, enzyme)
        new_fwd = rc_strand_cut(new_fwd, rev_site_index, enzyme)

        # create new reverse seq
        new_rev = self.rev
        new_rev = same_strand_cut(new_rev, rev_site_index, enzyme)
        new_rev = rc_strand_cut(new_rev, fwd_site_index, enzyme)

        return Fragment(new_fwd, new_rev)

    def get_overhangs(self):
        fwd_start_oh, fwd_end_oh = get_seq_overhangs(self.fwd)
        rev_start_oh, rev_end_oh = get_seq_overhangs(self.rev)
        return Overhang_set(fwd_start_oh, fwd_end_oh, rev_start_oh, rev_end_oh)

    def translate(self):
        dna = self.fwd.upper()
        start_index = re.search('ATG', dna).start()
        seq = Seq.Seq(dna[start_index:])
        return str(seq.translate())

    def __repr__(self):
        """ Temporary function for testing, works with bsaI specifically.
        """
        rev = self.rev[::-1]
        if self.fwd[:4].islower():
            return f'fwd: {self.fwd}\nrev:     {rev}'
        return f'fwd: {self.fwd}\nrev: {rev}'


def check_seq(frags):
    CDN_frags, AA_frags = zip(*frags)
    AA_seq = ''.join(AA_frags)
    building_frag = CDN_frags[0].digest(bsaI)
    for next_frag in CDN_frags[1:]:
        next_frag = next_frag.digest(bsaI)
        building_frag = concatenate_frags(building_frag, next_frag)
    CDN_seq = building_frag.translate()
    assert AA_seq == CDN_seq or (CDN_seq[0] == 'M' and AA_seq == CDN_seq[1:])
    assert not [m.start() for m in re.finditer('GGTCTC', building_frag.fwd)]
    assert not [m.start() for m in re.finditer('GGTCTC', building_frag.rev)]


def parent_frag_tuple(sr_id):
    """ Return a tuple (parent name, frag name) from the sr_id. Needs to
    handle case where parent name has underscores."""
    *name_list, frag = sr_id.split('_')
    name_list = [x + '_' for x in name_list[:-1]] + [name_list[-1]]
    parent_name = ''.join(name_list)
    return parent_name, frag


def verify_fragments(order_srs, AA_srs, lib_bps):
    order_frags = {parent_frag_tuple(sr.id): Fragment.from_SeqRecord(sr)
                   for sr in order_srs}
    order_frags = {(p_name, int(blk[-1])-1): v for (p_name, blk), v
                   in order_frags.items()}
    block_dict = {i: (lib_bps[i], b) for i, b in enumerate(lib_bps[1:])}
    AA_frags = {}
    for sr in AA_srs:
        seq = str(sr.seq)
        for b, (start, end) in block_dict.items():
            AA_frags.update({(sr.id, b): seq[start: end].replace('-', '')})

    combined_frags = {}
    for k, v in order_frags.items():
        combined_frags.update({k: (v, AA_frags[k])})

    pools = {}
    for (parent, block), v in combined_frags.items():
        if block in pools:
            pools[block].append((parent, v))
        else:
            pools[block] = [(parent, v)]
    pools = [pools[i] for i in pools]

    total_chimeras = len(pools[0]) ** len(pools)
    progress = general_tools.ProgressOutput(total_chimeras)
    for j, combo in enumerate(itertools.product(*pools)):
        progress.update(j)
        seq_frags = [i[1] for i in combo]
        check_seq(seq_frags)
    print('all sequences good')
