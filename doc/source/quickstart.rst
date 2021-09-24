
.. _quickstart:

**********
Quickstart
**********

0. Prequisites
--------------

To get started, you'll need a FASTA file with one or more parent amino acid sequences. If desired, we can find additional parents with BLAST.

Optionally, you can provide a PDB structure file, but otherwise we'll find one that matches the first parent sequence you provided. See :class:`schemarecomb.PDBStructure` for more details.

.. note::

    This guide assumes you're using Python 3.9 on Linux. If you use MacOS, things will probably work the same, but no guarantees. If you have Windows, I recommend you use the `Windows Subsystem for Linux  <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_, but again, no guarantees. Please raise an issue on the schemarecomb GitHub page if you have OS difficulty.
   

1. Install schemarecomb
-----------------------

schemarecomb is available on pip::

    $ pip install schemarecomb

Or you can install schemarecomb from source. 

In a Python script, import schemarecomb::

    import schemarecomb


2. Make a ParentSequences
-------------------------

Load your parent FASTA file and find additional parents if needed. For this example, we'll use beta-glucosidase (bgl3, PDB ID 1GNX)::

    parent_fn = 'bgl3.fasta'
    p_aln = schemarecomb.ParentSequences.from_fasta(parent_fn)
    p_aln.obtain_seqs(num_final_seqs=4, desired_identity=0.7)
    p_aln.get_PDB()

After running, p_aln is a ParentSequences with four parents that have about 70% pairwise identity and the closest PDB structure.

See :class:`schemarecomb.ParentSequences` for more options. Viewing :class:`schemarecomb.PDBStructure` may also be helpful.


3. Run the SCHEMA-RASPP algorithm
---------------------------------

SCHEMA-RASPP finds potential libraries and calculates the probability of Golden Gate assembly for each::

    libraries = schemarecomb.generate_libraries(p_aln, 6)

This finds libraries with six blocks (five breakpoints).

See :func:`schemarecomb.generate_libraries` for more options.


4. Select and save a library
----------------------------

Select the library with highest mutation_rate - energy and save generated DNA blocks::

    best_lib = max(libraries, key=lambda x: x.mutation_rate - x.energy)
    SeqIO.write(best_lib.dna_blocks, 'bgl3_library_dna.fasta', 'fasta')

The DNA fragments in FASTA format in a file named "bgl3_library_dna.fasta". These fragments are ready to order and assemble with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. You can simulate the Golden Gate reaction using SnapGene.

See :class:`schemarecomb.Library` for more options.
