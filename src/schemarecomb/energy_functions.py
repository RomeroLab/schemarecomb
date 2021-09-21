# energy_functions.py

"""Chimeric library energy calculations.

Library energy should be correlated with the functional fraction of proteins in
the library. Thus, a library with high energy is likely to have a smaller
number of active chimeras than a low energy library.

Currently, only the SCHEMA energy function is implemented for SCHEMA-RASPP
calculation (Endelman et al. 2004). This module is subject to change in future
versions.

"""

from abc import abstractmethod
from abc import ABC
from decimal import Decimal
from importlib import import_module
from itertools import product

import numpy as np

import schemarecomb


class EnergyFunction(ABC):
    """Abstract class for making energy functions."""
    # TODO: Implement general energy function. (v0.2.0)

    def __init__(self, parents: schemarecomb.ParentSequences):
        self.parents = parents

    @abstractmethod
    def block_energy(self, start: int, end: int) -> Decimal:
        raise NotImplementedError

    @abstractmethod
    def increment_block_energy(self, start: int, new_end: int) -> Decimal:
        raise NotImplementedError

    @property
    @abstractmethod
    def import_mod_cls(self) -> tuple[str, str]:
        raise NotImplementedError


def build_from_str(mod: str, name: str, parents: schemarecomb.ParentSequences):
    """Build EnergyFunction subclass instance from name and ParentSequences."""
    e_function_type = getattr(import_module(mod), name)
    return e_function_type(parents)


class SCHEMA(EnergyFunction):
    r"""Calculates the energy (approximate fraction functional) of a library.

    Given a PDB structure and a multiple sequence alignment, the SCHEMA energy
    of a recombinant sequence is the number of pairwise interactions broken by
    recombination. Expressed as an equation (Endelman et al. 2004 eq. 7):

    .. math::

        S(s) = \sum_i \sum_{j>i} C(i, j) * D(s[i], s[j])

    where:
        :math:`s` is the recombinant sequence.

        :math:`C(i, j) = 1` if residues i and j are within 4.5 Angstroms in the
        protein structure, else 0.

        :math:`D(s[i], s[j]) = 0` if the amino acids at positions i and j in
        the recombinant sequence are also present together in any parent,
        else 1.

    For example, if "AEH" and "GQH" are aligned parents with contacts (0, 1)
    and (1, 2), chimera "AQH" has a SCHEMA energy of 1.0. Residues 'A' and 'Q'
    are not found together in any parents at positions 0 and 1, respectively,
    and there is a contact at these positions, so
    :math:`C(0, 1) * D(s[0], s[1]) = 1`. The 'Q' and 'H', however, do not
    contribute to the SCHEMA energy. Although :math:`C(1, 2) = 1` because
    residues 1 and 2 are contacting, 'Q' and 'H' are found together at these
    positions in the second parent, so :math:`D(s[1], s[2]) = 0`.

    See Voigt et al. 2002 and Silberg et al. 2004 for more information about
    SCHEMA energy.

    The E_matrix attribute is a square matrix with the same size as the length
    of alignment. :math:`E_{rj}` is Endelman et al. 2004 eq. 6. with the r
    and t sums removed and evaluated for average SCHEMA energy for a
    specific r and t:

    .. math::

        E_{rj} = \frac{1}{p^2} \sum_{p \in parents} \sum_{q \in parents}
        C(r, j) * D(p[r], q[j]).

    Then eq. 6 can be evaluated for a new block starting at position X_k by
    summing over :math:`E` as in the block_energy method:

    .. math::

        \langle E \rangle_{(X_1, X_2, ..., X_{k-1}, X_k)}
        \langle E \rangle_{(X_1, X_2, ..., X_{k-1})} = \sum_{r=X_{k-1}}^{X_k-1}
        \sum_{t=X_k}^{N-1} E_{rj}.

    Note the sums are decremented by one compared to eq. 6 as it appears in
    Endelman et al 2004. This was done to better reflect Python indexing.

    Finally, a computational speed up is to sequentially calculate the energy
    of a new block starting at each :math:`X_k` given the energy of a new block
    at position :math:`X_k - 1`. This is implemented in the
    increment_block_energy method.

    Parameters:
        parents: Parent sequences for recombination. Must be aligned and have
            pdb_structure attribute. Note that the SCHEMA object is initialized
            for the ParentSequences at that point in time, so if parents
            changes, the SCHEMA object will NOT change also.

    Attributes:
        E_matrix (np.ndarray): Energy matrix corresponding to input
            ParentSequences' alignment and contacts.
        parents: Parent sequences for recombination.

    Raises:
        ValueError: If input ParentSequences is not aligned or does not have
            the pdb_structure attribute on initialization.

    """
    def __init__(self, parents: schemarecomb.ParentSequences):
        super().__init__(parents)
        try:
            alignment: list[tuple[str, ...]] = parents.alignment
        except AttributeError:
            raise ValueError('Input ParentSequences must be aligned.')

        try:
            pdb_structure = parents.pdb_structure
        except AttributeError:
            raise ValueError('Input ParentSequences must have a '
                             'pdb_structure.')

        E_matrix = np.zeros((len(alignment), len(alignment)))

        # Only need to evaluate at contacts because otherwise C_{i,j} = 0.
        for i, j in pdb_structure.contacts:
            # Get the ordered parent amino acids at i and j.
            AAs_i, AAs_j = alignment[i], alignment[j]

            # Parent amino acid pair might be redundant.
            parental = set(zip(AAs_i, AAs_j))

            # For each possible recombinant calculate D_ij. E_matrix entry
            # is this sum over the number of possible recombinants.
            broken = sum(1 for combo in product(AAs_i, AAs_j) if combo
                         not in parental)
            E_matrix[i, j] = broken / (len(AAs_i)**2)

        self.E_matrix = E_matrix
        # self.parents = parents

    def block_energy(self, start: int, end: int) -> Decimal:
        r"""Average SCHEMA energy of adding a block to a library.

        This is Endelman et al. 2004 eq 6. Move the r and t sums to the outside
        of the equation. Then the inner part of the equation is E_matrix,
        described in the class overview. The summed sliced returned by this
        function thus corresponds to evaluating the outside r and t sums.

        Parameters:
            start: Index of the first residue the block that was previously the
                last block. This corresponds to :math:`x_k-1` in the class
                overview.
            end: Index of the first residue in the new block. Last residue in
                the previous block is at end-1. end corresponds to
                :math:`x_k` in the class overview.

        Returns:
            Average SCHEMA energy of adding the block.
        """
        energy = self.E_matrix[start:end, end:].sum()
        return Decimal(energy)

    def increment_block_energy(self, start: int, new_end: int) -> Decimal:
        """Difference in average SCHEMA energy when block size is decremented.

        This is equivalent to
            block_energy(start, new_end) - block_energy(start, new_end-1).

        This function is faster than separate block_energy calculations.
        Therefore, you can calculate the energy for the first node in the
        graph using block_energy, then calculate the subsequent node energies
        using this function.

        Parameters:
            start: Index of the first residue in the block before the new
                block. This corresponds to :math:`x_k-1` in the class overview.
            new_end: Index of the first residue in the new block, after the
                new block size is decremented. This corresponds to
                :math:`x_k` in the class overview.

        Returns:
            Difference in average SCHEMA energy between the new block starting
            at new_end and starting at new_end-1.
        """
        old_end = new_end - 1
        neg = self.E_matrix[start:old_end, old_end].sum()
        pos = self.E_matrix[old_end, new_end:].sum()
        return Decimal(pos - neg)

    @property
    def import_mod_cls(self) -> tuple[str, str]:
        """Save module and class so we know how to construct it later."""
        return self.__module__, type(self).__name__
