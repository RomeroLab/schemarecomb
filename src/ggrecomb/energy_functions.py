from itertools import product

import numpy as np

from ggrecomb.parent_alignment import ParentAlignment


class EnergyFunction:
    pass


class SCHEMA(EnergyFunction):
    """SCHEMA energy.

    Given a PDB structure and a multiple sequence alignment, the SCHEMA energy
    of a recombinant sequence is the number of pairwise interactions broken by
    recombination. Expressed as an equation (Endelman et al. 2004 eq. 7):

        E = sum_i sum_{j>i} C_{i,j} D_{i,j}

        where:
            C_{i,j} = 1 if residues i and j are within 4.5 Angstroms, else 0.
            D_{i,j} = 0 if the amino acids at positions i and j in the
                recombinant sequence are also present together in any parent,
                else 0.

    Example: If "AEH" and "GQH" are aligned parents with contacts (0, 1) and
        (1, 2), chimera "AQH" has a SCHEMA energy of 1.0 because residues 'A'
        and 'Q' are not found in any parents at positions 0 and 1,
        respectively, and there is a contact at these positions. Although
        residues 1 and 2 are contacting, 'Q' and 'H' are found together at
        these positions in the second parent, so D_{0,1} = 0 and doesn't
        contribute to the SCHEMA energy.

    See Voigt et al. 2002 and Silberg et al. 2004 for more information.

    Public methods:
        block_energy: Average SCHEMA energy of a section of a parent multiple
            sequence alignment.
        increment_block_energy: Average SCHEMA energy difference between a
            block and a block decreased by one residue.
    """
    def __init__(self, pa: ParentAlignment):
        """Generate the SCHEMA energy matrix for block calculation.

        E_matrix is a square matrix with the same size as the length of
        alignment. E_matrix_{r,j} is Endelman et al. 2004 eq. 6. with the r and
        t sums removed and evaluated for average SCHEMA energy for a specific r
        and t. Then eq. 6 can be evaluated for a given block by summing over
        E_matrix_{r,t}.

        Args:
            pa: Parent sequence alignment. Must have pdb_structure attribute.
        """
        alignment = list(zip(*pa.aligned_sequences))
        self.E_matrix = np.zeros((len(alignment), len(alignment)))

        # Only need to evaluate at contacts because otherwise C_{i,j} = 0.
        for i, j in pa.pdb_structure.contacts:
            # Get the ordered parent amino acids at i and j.
            AAs_i, AAs_j = alignment[i], alignment[j]

            # Parent amino acid pair might be redundant.
            parental = set(zip(AAs_i, AAs_j))

            # For each possible recombinant calculate D_ij. E_matrix entry
            # is this sum over the number of possible recombinants.
            broken = sum(1 for combo in product(AAs_i, AAs_j) if combo
                         not in parental)
            self.E_matrix[i, j] = broken / (len(AAs_i)**2)

    def block_energy(self, start: int, end: int) -> float:
        """Calculate the average SCHEMA energy of a block.

        This is Endelman et al. 2004 eq 6.

        Args:
            start: Index of the first residue in the block.
            end: Index of the first residue in the next block. Last residue in
                the current block is at end-1.

        Returns:
            Average SCHEMA energy of block [start, end].
        """
        energy = self.E_matrix[start:end, end:].sum()
        return energy

    def increment_block_energy(self, start: int, new_end: int) -> float:
        """Difference in average SCHEMA energy when block size is incremented.

        This is equivalent to
            block_energy(start, new_end) - block_energy(start, new_end-1).

        This function is faster than separate block_energy calculations.
        Therefore, you can calculate the energy for the first node in the
        graph using block_energy, then calculate the subsequent node energies
        using this function.

        Args:
            start: Index of the first residue in the block.
            end: Index of the first residue in the next block, after the
                current block size is incremented.

        Returns:
            Difference in average SCHEMA energy between block [start, new_end]
                and block [start, new_end-1].
        """
        old_end = new_end - 1
        neg = self.E_matrix[start:old_end, old_end].sum()
        pos = self.E_matrix[old_end, new_end:].sum()
        return pos - neg
