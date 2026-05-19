"""Functions and classes that handle :class:`~schemarecomb.Library` instances.

This module provides the definition of :class:`schemarecomb.Library`, which
represents a library of chimeric proteins, and functions for analyzing and
manipulating libraries.

"""

from collections import deque
from collections.abc import MutableMapping
from dataclasses import dataclass
from decimal import Decimal
from functools import cached_property
from importlib import metadata
from itertools import combinations
import json
from re import fullmatch
from random import choice
from typing import Optional
from typing import Union

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import schemarecomb
from schemarecomb.breakpoints import BreakPoint
from schemarecomb.breakpoints import block_indices
from schemarecomb.breakpoints import Overhang
from schemarecomb import energy_functions
from schemarecomb.energy_functions import EnergyFunction
from schemarecomb.restriction_enzymes import RestrictionEnzyme


# List of 31 common codons {+stop and gap) in E. coli.
_AA_C31: dict[str, set[str]]
_AA_C31 = {'A': {'GCT', 'GCA'}, 'R': {'CGT', 'CGA'}, 'N': {'AAT'},
           'D': {'GAT'}, 'C': {'TGT'}, 'Q': {'CAA', 'CAG'}, 'E': {'GAA'},
           'G': {'GGT'}, 'H': {'CAT', 'CAC'}, 'I': {'ATT', 'ATC'},
           'L': {'TTA', 'TTG', 'CTA'}, 'K': {'AAA'}, 'M': {'ATG'},
           'F': {'TTT'}, 'P': {'CCT', 'CCA'}, 'S': {'AGT', 'TCA'},
           'T': {'ACA', 'ACT'}, 'W': {'TGG'}, 'Y': {'TAT'},
           'V': {'GTT', 'GTA'}, '*': {'TGA'}, '-': {'---'}}


class MutationRateCache(MutableMapping):
    """Caches calculated mutation rates for future mutation calculations.

    For certain breakpoints, shifting it by one residue will not change the
    average mutation rate of otherwise identical libraries. Therefore, we can
    cache every rate that is calculated, such that we can use it for future
    libraries that are known to differ only in these redundant breakpoints.
    This is necessary because the average mutation rate calculation takes a
    significant amount of time and is otherwise difficult to speed up.

    You can treat this class as a dictionary with keys of breakpoint position
    tuples and values of mutation_rates. The breakpoint positions will be
    automatically converted to breakpoint grouping index. For example, suppose
    otherwise identical libraries with breakpoints at positions 2 or 3 are
    known to have the same mutation rate. Then::

        >>> from schemarecomb.libraries import MutationRateCache
        >>> bp_to_group = {1: 1, 2: 2, 3: 2, 4: 3}  # 2 and 3 in the same group
        >>> mr_cache = MutationRateCache(bp_to_group)
        >>> mr_cache[(1, 2, 4)] = 23.4  # avg mutation rate of lib (1, 2, 4)
        >>> # We know the (1, 3, 4) avg mutation rate without calculation.
        >>> assert mr_cache[(1, 3, 4)] == 23.4

    Parameters:
        bp_to_group_bp: Mapping from breakpoint position to the index of group
            of redundant breakpoints that contain it.

    """

    def __init__(self, bp_to_group_bp: dict[int, int]):
        self._bp_to_group_bp = bp_to_group_bp
        self._group_to_m: dict[tuple[int, ...], float] = {}

    @classmethod
    def from_parents(cls, parents: schemarecomb.ParentSequences,
                     valid_bps: list[BreakPoint]) -> 'MutationRateCache':
        """Construct from ParentSequences and a list of valid breakpoints.

        Parameters:
            parents: Parent alignment used to calculate breakpoint groupings.
            valid_bps: Valid breakpoints to group together.

        """
        # Group the redundant breakpoints.
        breakpoint_groups = []

        # Iterate over the internal (not added start or end) breakpoints.
        bps = iter(sorted([bp.position for bp in valid_bps]))
        curr_bps = [next(bps)]
        for bp in bps:
            # If previous bp is not adjacent to current or has more than one
            # amino acid, current belongs to a new group.
            # TODO: Test new group logic.
            if bp != curr_bps[-1] + 1 or \
               len(set(parents.alignment[bp-1])) != 1:
                breakpoint_groups.append(curr_bps)
                curr_bps = [bp]
            else:
                curr_bps.append(bp)
        breakpoint_groups.append(curr_bps)

        # Add 0 and len(parents.alignment) if not already included, for safety.
        # Nothing gets hurt if these are never accessed. Libraries that don't
        # have overhangs at 0 and aln_len may or may not have these indices in
        # their representations.
        bp_positions = {bp.position for bp in valid_bps}
        if 0 not in bp_positions:
            breakpoint_groups.append([0])
        aln_len = len(parents.alignment)
        if aln_len not in bp_positions:
            breakpoint_groups.append([aln_len])

        # Get the mapping from each breakpoint index to the breakpoint's group
        # index.
        bp_to_group = {}
        for group_index, bp_group in enumerate(breakpoint_groups):
            for bp in bp_group:
                bp_to_group[bp] = group_index

        return cls(bp_to_group)

    def __getitem__(self, key: tuple[int, ...]) -> float:
        # convert bps in key to group bps
        # if group bps in group_to_m, return
        # else return None, then cache later
        group_bps = tuple(self._bp_to_group_bp[bp] for bp in sorted(key))

        # Check that all elements of group_bos are unique.
        if len(set(group_bps)) != len(group_bps):
            raise KeyError('At least two input breakpoints are redundant.')

        try:
            return self._group_to_m[group_bps]
        except KeyError:
            raise KeyError('Non-redundant breakpoint set not in cache.')

    def __setitem__(self, key: tuple[int, ...], value: float) -> None:
        # convert bps in key to group_bps
        # set self.group_to_m[group_bps] = value
        group_bps = tuple(self._bp_to_group_bp[bp] for bp in sorted(key))

        # Check that all elements of group_bos are unique.
        if len(set(group_bps)) != len(group_bps):
            raise KeyError('At least two input breakpoints are redundant.')

        self._group_to_m[group_bps] = value

    def __delitem__(self, key):
        raise TypeError('Deletion is not supported for MutationRateCache.')

    def __iter__(self):
        return iter(self._group_to_m)

    def __len__(self):
        return len(self._group_to_m)


def _sequence_mutations(seq1, seq2):
    """Hamming distance between seq1 and seq2."""
    assert len(seq1) == len(seq2)
    return sum(1 for a1, a2 in zip(seq1, seq2) if a1 != a2)


def average_mutations(
    breakpoints: list[BreakPoint],
    parents: schemarecomb.ParentSequences,
    mr_cache: Optional[MutationRateCache] = None
) -> Decimal:
    """Calculates the average number of mutations in a chimeric library.

    Parameters:
        breakpoints: Library recombination points.
        parents: Parent alignment used to calculate mutation rate.
        mr_cache: Caches mutation rates to speed up calculations.

    Returns:
        Average number of mutations from chimera to nearest parent for each
            chimera in library.

    """

    bp_indices = tuple(sorted(bp.position for bp in breakpoints))
    if mr_cache is not None and \
       (mut_rate := mr_cache.get(bp_indices, None)) is not None:
        return mut_rate

    seqs = list(zip(*parents.alignment))

    # Get mutations of each block.
    blmuts = []
    aln_len = len(parents.alignment)
    block_inds = block_indices(breakpoints, aln_len)
    for blk_start, blk_end in block_inds:
        muts = {i: {i: 0} for i, _ in enumerate(seqs)}
        for (i, s1), (j, s2) in combinations(enumerate(seqs), 2):
            s1_slice = s1[blk_start:blk_end]
            s2_slice = s2[blk_start:blk_end]
            m = _sequence_mutations(s1_slice, s2_slice)
            muts[i][j] = m
            muts[j][i] = m
        blmuts.append(muts)

    # Calculate parameters for mutation calculation.
    num_blocks = len(block_inds)
    num_parents = len(seqs)

    # Perform depth first traversal over blocks to calculate M_sum. This
    # is faster than itertools.product because it avoids redundant
    # summation.
    stack = deque((1, [blmuts[0][p1][p2] for p2 in range(num_parents)])
                  for p1 in range(num_parents))
    mut_rate_sum = 0
    while stack:
        i, par_muts = stack.pop()

        par_iter = list(zip(par_muts, range(num_parents)))
        block_muts = blmuts[i]
        for p1 in range(num_parents):
            p_muts = block_muts[p1]
            new_part_muts = [m + p_muts[p2] for m, p2 in par_iter]
            if i + 1 == num_blocks:
                mut_rate_sum += min(new_part_muts)
            else:
                stack.append((i + 1, new_part_muts))

    mut_rate = Decimal(mut_rate_sum) / Decimal(num_parents**num_blocks)

    if mr_cache is not None:
        mr_cache[bp_indices] = mut_rate

    return mut_rate


class LibraryConfig:
    """Holds parameters shared across a collection of Library instances.

    Passed into :func:`schemarecomb.Library.calc_from_config` constructor.

    Parameters:
        energy_function: The EnergyFunction used the calculate to energy
            attribute of a Library. Note that this attribute is purely for
            recordkeeping; no computation is done with energy_function in
            Library.calc_from_config. Therefore, it is up to the user to ensure
            that energy_function correctly reflects the EnergyFunction used
            calculate energy.
        gg_enzyme: The enzyme used to calculate gg_prob and optimize over
            possible gg_overhangs. This should correspond to the restriction
            enzyme to be used to assemble the library.
        gg_threshold: Minimum gg_prob for a set of overhangs to be considered
            optimal. May be set < 1.0 for faster runtime, but then a Library's
            gg_overhangs is not guaranteed to be the set with maximum gg_prob.
        mr_cache: Optional cache used to speed up mutation rate calculations
            by reusing values from adjacent libraries.
        amino_to_cdn: Mapping of amino acids to codons. Forms a simple codon
            optimization.

    """

    def __init__(
        self,
        energy_function: EnergyFunction,
        gg_enzyme: RestrictionEnzyme,
        gg_threshold: Union[float, Decimal] = 1.0,
        mr_cache: Optional[MutationRateCache] = None,
        amino_to_cdn: Optional[dict[str, set[str]]] = None,
    ):
        self.energy_function = energy_function
        self.gg_enzyme = gg_enzyme
        self.gg_threshold = gg_threshold
        self.mr_cache = mr_cache
        if amino_to_cdn is None:
            self.amino_to_cdn = _AA_C31


@dataclass
class _Library:
    """A chimeric protein library.

    This class is not usually instantiated by users, but instead by another
    package function. If generating many libraries, the simplest construction
    method is via the :meth:`schemarecomb.calc_from_config` constructor, which
    uses a :class:`~schemarecomb.libraries.LibraryConfig` to consolidate or
    compute all parameters except breakpoints and energy.

    A library consists of a parent alignment and a series of breakpoints, which
    are the indices of the alignment where the parent sequences are recombined
    with Golden Gate assembly to form chimeras. A "block" is defined as the
    sequence between adjacent breakpoints, such that each chimeric sequence is
    composed of sequential blocks each taken from a parent. For example, a
    library with 2 breakpoints and 3 parents has a chimera that may be
    identified as "201": the first block from parent 2, the second block from
    parent 0, and the last block from parent 1.

    Each library is assigned metrics of energy and mutation_rate that estimate
    the fraction of functional recombinants and sequence diversity,
    respectively. The exact metrics depend on the method used to optimize over
    the recombinant space. For example, the SCHEMA-RASPP algorithm averages the
    SCHEMA energy and minimum number of mutations from a parent sequence over
    every recombinant protein in the library. See Endelman et al. 2004 for more
    details about SCHEMA-RASPP.

    Each library is also assigned a "Golden Gate probability" that estimates
    the efficiency of library assembly given the optimal set of Golden Gate
    overhangs chosen for the given library. This value is computed (and
    overhangs chosen) based on data from Gregory Lohman and collaborators at
    NEB (Potapov et al. 2018). See Pryor, Potapov, et al. 2020 for details on
    the biology and mathematics behind this technique. (Note that the
    calculation was developed independently by the authors of schemarecomb in
    2019 after viewing a talk by Dr. Lohman.)

    See :mod:`schemarecomb.libraries` for classes and functions that handle
    Library instances.

    The attributes for this class are the same as the parameters, with the
    addition of max_block_len and min_block_len attributes that give the
    maximum and minimum block length in the library, respectively.

    You might want to add additional DNA bases to the ends of the fragments, as
    this may improve restriction enzyme efficiency.

    Note:
        This class may produce unintended Golden Gate sites when generating
        DNA. The current version does not check for this. Simulate restriction
        enzyme cutting and Golden Gate Assembly in a program like Benchling or
        Snapgene, then change codons as necessary.

    Parameters:
        breakpoints: Alignment indices where the parent sequences are
            recombined to form the recombinants in the Library.
        energy: Estimation of the fraction of functional recombinants relative
            to other libraries from the same parent alignment. High energy
            libraries are likely to have a larger proportion of active enzymes.
            The interpretation of this value depends on the EnergyFunction
            used to calculate it, which must be passed in as the
            energy_function parameter.
        energy_function: The EnergyFunction used the calculate the energy
            attribute.
        mutation_rate: Average mutational distance to closest parent for each
            recombinant in the library.
        gg_prob: Golden Gate efficiency using gg_overhangs, calculated with
            gg_enzyme.
        gg_overhangs: Collection of overhangs found with (near) maximum
            gg_prob. Element indices correspond to the same index in
            breakpoints, e.g. gg_overhangs[i] is the overhang at
            breakpoints[i].position. Note that there may exist a valid overhang
            set for this library with better gg_prob depending on the heuristic
            used to optimize gg_prob.
        gg_enzyme: The enzyme used to calculate gg_prob and optimize over
            possible gg_overhangs. This should correspond to the restriction
            enzyme to be used to assemble the library.
        amino_to_cdn: Mapping of amino acids to codons. Forms a simple codon
            optimization.

    Attributes:
        All parameters for this class are also attributes.
        block_indices (list[tuple[int, int]]): The start and end indicies of
            each block. See :func:`~schemarecomb.breakpoints.block_indices`
            for more information.
        max_block_len (int): Size of the largest block in the library.
        min_block_len (int): Size of the smallest block in the library.
        dna_blocks (list[list[str]]): DNA sequences for the chimeric blocks,
            with Golden Gates sites at the ends. Each inner list corresponds to
            one of the parents in energy_function.parents. Each element in an
            inner list is a DNA sequence that translates to amino acid blocks
            from the parent.
        dna_blocks (list[SeqRecord]): Chimeric DNA blocks with id
            '<parent_name>_block-<block number>'. Groups of chimeric blocks
            with compatible block numbers (from 0 to number of blocks,
            inclusive and all unique) may be Golden Gate Assembled to form DNA
            that may be transcribed and translated into chimeric proteins.

    """

    breakpoints: list[BreakPoint]
    energy: Decimal
    energy_function: EnergyFunction
    mutation_rate: Decimal
    gg_prob: Decimal
    gg_overhangs: list[Overhang]
    gg_enzyme: RestrictionEnzyme
    amino_to_cdn: dict[str, set[str]]

    @cached_property
    def block_indices(self):
        aln_len = len(self.energy_function.parents.alignment)
        return block_indices(self.breakpoints, aln_len)

    @cached_property
    def min_block_len(self):
        return min(end - start for start, end in self.block_indices)

    @cached_property
    def max_block_len(self):
        return max(end - start for start, end in self.block_indices)

    @property
    def _schemarecomb_version(self):
        return metadata.version('schemarecomb')

    @property
    def dna_blocks(self):
        """
        Generates DNA blocks for Golden Gate assembly.

        Finds a single 'design' codon/amino acid combination chosen by
        the optimizer for each junction by testing all possibilities for all
        homologues at the correct sequence position, ensuring consistent
        overhangs.
        """
        try:
            return self._dna_blocks
        except AttributeError:
            parents = self.energy_function.parents
            breakpoints = self.breakpoints
            blk_indices = self.block_indices
            overhangs = self.gg_overhangs.copy()

            breakpoint_map = breakpoints.copy()
            if blk_indices[0] != breakpoints[0].position:
                overhangs.insert(0, None)
                breakpoint_map.insert(0, None)

            if blk_indices[-1] != breakpoints[-1].position:
                overhangs.append(None)
                breakpoint_map.append(None)

            start_end_overhangs = list(zip(overhangs, overhangs[1:]))
            start_end_breakpoints = list(
                zip(breakpoint_map, breakpoint_map[1:]))

            resenz = self.gg_enzyme
            oh_len = resenz.overhang_len
            dna_blocks = []

            for parent_sr in parents.records:
                parent_aa = parent_sr.seq
                par_blocks = []
                for blk_num, positions in enumerate(blk_indices):
                    oh_start, oh_end = start_end_overhangs[blk_num]
                    bp_start, bp_end = start_end_breakpoints[blk_num]
                    start, end = positions

                    aa_block = [aa for aa
                                in parent_aa[start:end] if aa != '-']
                    aa_start_parent, *aa_block_internal, \
                        aa_end_parent = aa_block
                    cdn_options = [self.amino_to_cdn[aa] for aa in
                                   aa_block_internal]
                    dna_block_internal = ''.join(
                        choice(tuple(cdn_set)) for cdn_set in cdn_options)

                    # Handle start of the block
                    if oh_start is None:
                        beg_dna = choice(
                            tuple(self.amino_to_cdn[aa_start_parent]))
                    else:
                        end_dot = 6 - oh_start.ind - oh_len
                        pattern = (oh_start.seq + '.' * end_dot)[-3:]

                        pos = bp_start.position
                        possible_aas = {p.seq[pos] for p in parents.records if
                                        p.seq[pos] != '-'}

                        matching_codon = ''
                        for aa in possible_aas:
                            if aa in self.amino_to_cdn:
                                for codon in self.amino_to_cdn[aa]:
                                    if fullmatch(pattern, codon):
                                        matching_codon = codon
                                        break
                            if matching_codon:
                                break

                        if not matching_codon:
                            # Should never happen
                            raise ValueError(
                                f"Could not find a matching codon for pattern "
                                f"'{pattern}' at position {pos}."
                                f" Searched AAs: {possible_aas}")

                        if end_dot:
                            overhang_amino = oh_start.seq.lower()\
                                             + matching_codon[-end_dot:]
                        else:
                            overhang_amino = oh_start.seq.lower()
                        beg_dna = resenz.for_restriction_site + overhang_amino

                    # Handle end of the block
                    if oh_end is None:
                        end_dna = choice(
                            tuple(self.amino_to_cdn[aa_end_parent]))
                    else:
                        begin_dot = oh_end.ind
                        pattern = ('.' * begin_dot + oh_end.seq)[:3]

                        #: The amino acid being modified is at position-1 of
                        # the breakpoint.
                        pos = bp_end.position - 1

                        possible_aas = {p.seq[pos] for p in parents.records if
                                        p.seq[pos] != '-'}

                        matching_codon = ''
                        for aa in possible_aas:
                            if aa in self.amino_to_cdn:
                                for codon in self.amino_to_cdn[aa]:
                                    if fullmatch(pattern, codon):
                                        matching_codon = codon
                                        break
                            if matching_codon:
                                break

                        if not matching_codon:
                            # Should never raise
                            raise ValueError(
                                f"Could not find a matching codon for pattern"
                                f" '{pattern}' at position {pos + 1} "
                                f"(AA at pos {pos}). "
                                f"Searched AAs: {possible_aas}")

                        if begin_dot:
                            amino_overhang = matching_codon[
                                             :begin_dot] + oh_end.seq.lower()
                        else:
                            amino_overhang = oh_end.seq.lower()
                        end_dna = amino_overhang + resenz.rev_restriction_site

                    block_seq = Seq(beg_dna + dna_block_internal + end_dna)
                    block_name = f'{parent_sr.id}_block-{blk_num}'
                    block_sr = SeqRecord(block_seq, id=block_name,
                                         name=block_name,
                                         description=f'{parent_sr.description}'
                                                     f' | block {blk_num}')
                    par_blocks.append(block_sr)

                dna_blocks.extend(par_blocks)

            self._dna_blocks = dna_blocks
            return dna_blocks

    def find_best_overhangs(self) -> None:
        """Find the overhangs with the true maximum gg_prob.

        Resets the gg_prob and gg_overhangs attributes. This method can be used
        if a gg_prob threshold <1.0 was used during construction to speed up
        calculation. After this library is selected, call this method to find
        the truly optimal set of gg_overhangs, so that the (possibly) improved
        overhangs are used in the generated DNA sequences.

        """
        # No thresh parameter => thresh=1.0 => optimal overhangs.
        best_prob, best_ohs = self.gg_enzyme.max_gg_prob(self.breakpoints)
        self.gg_prob = best_prob
        self.gg_overhangs = best_ohs

    @classmethod
    def calc_from_config(
        cls,
        breakpoints: list[BreakPoint],
        energy: Decimal,
        config: LibraryConfig
    ) -> '_Library':
        """Calculate average mutation rate and gg_prob during construction.

        Useful for when constructing many libraries from the same
        :class:`~schemarecomb.ParentSequences`. Energy must be precalculated.

        Parameters:
            breakpoints: The breakpoints that define the new library.
            energy: The energy of the new library.
            config: Configuration shared between all other libraries.
        """
        parents = config.energy_function.parents

        mut_rate = average_mutations(breakpoints, parents, config.mr_cache)

        gg_thresh = config.gg_threshold
        if isinstance(gg_thresh, float):
            gg_thresh = Decimal(gg_thresh)

        gg_prob, gg_overhangs = config.gg_enzyme.max_gg_prob(
            breakpoints,
            gg_thresh
        )

        return cls(
            breakpoints,
            energy,
            config.energy_function,
            mut_rate,
            gg_prob,
            gg_overhangs,
            config.gg_enzyme,
            config.amino_to_cdn
        )

    def __setattr__(self, name, value) -> None:
        super().__setattr__(name, value)

        # If gg_overhangs changes, cached _dna_blocks are no longer valid.
        if name == 'gg_overhangs' and hasattr(self, '_dna_blocks'):
            del self._dna_blocks

    def to_json(self) -> str:
        """Convert instance to a JSON-formatted string.

        Return:
            Instance converted to a JSON string.

        """
        bp_list = [[bp.position, [[oh.ind, oh.seq] for oh in bp.overhangs]]
                   for bp in self.breakpoints]

        amino_to_cdn = {aa: tuple(cdns) for aa, cdns
                        in self.amino_to_cdn.items()}

        out_list = [
            bp_list,
            str(self.energy),
            self.energy_function.import_mod_cls,
            self.energy_function.parents.to_json(),
            str(self.mutation_rate),
            str(self.gg_prob),
            self.gg_overhangs,
            self.gg_enzyme.to_json(),
            amino_to_cdn,
            f'schemarecomb version: {self._schemarecomb_version}'
        ]
        return json.dumps(out_list)

    @classmethod
    def from_json(cls, in_json: str) -> '_Library':
        """Construct instance from JSON.

        Parameters:
            in_json: JSON-formatted string representing a Library.

        Return:
            ParentSequences instance created from in_json.

        """

        breakpoints_json, energy, *other_attrs = json.loads(in_json)
        e_function_tup, pa_json, *other_attrs = other_attrs
        mr, gg_prob, *other_attrs = other_attrs
        gg_overhangs, gg_enzyme_str, *other_attrs = other_attrs
        amino_to_cdn_tups, _ = other_attrs

        breakpoints = [BreakPoint(pos, [Overhang(*oh) for oh in ohs])
                       for pos, ohs in breakpoints_json]

        energy = Decimal(energy)

        parents = schemarecomb.ParentSequences.from_json(pa_json)

        ef_mod, ef_name = e_function_tup
        e_function = energy_functions.build_from_str(ef_mod, ef_name, parents)

        mr = Decimal(mr)
        gg_prob = Decimal(gg_prob)

        gg_overhangs = [Overhang(*gg_oh_tup) for gg_oh_tup in gg_overhangs]

        gg_enzyme = RestrictionEnzyme.from_json(gg_enzyme_str)

        amino_to_cdn = {aa: set(cdns) for aa, cdns
                        in amino_to_cdn_tups.items()}

        return cls(breakpoints, energy, e_function, mr, gg_prob, gg_overhangs,
                   gg_enzyme, amino_to_cdn)

    # TODO: Might be needed in raspp?
    '''
    def expand(self, e_diff: float, new_bp: dict[int, list[tuple[int, str]]]):
        new_e = self.energy + e_diff
        new_bps = self.breakpoints | new_bp
        return type(self)(new_e, new_bps, self.parent_alignment)
    '''
