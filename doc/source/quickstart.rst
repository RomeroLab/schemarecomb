
.. _quickstart:

**********
Quickstart
**********

0. Prequisites
--------------

To get started, you'll need a FASTA file with one or more parent amino acid sequences. If desired, we can find additional parents with BLAST.

Optionally, you can provide a PDB structure file, but otherwise we'll find one that matches the first parent sequence you provided. See :class:`ggrecomb.PDBStructure` for more details.

.. note::

    This guide assumes you're using Python 3.9 on Linux. If you use MacOS, things will probably work the same, but no guarantees. If you have Windows, I recommend you use the `Windows Subsystem for Linux  <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_, but again, no guarantees. Please raise an issue on the ggrecomb GitHub page if you have OS difficulty.
   

1. Install ggrecomb
-------------------

ggrecomb is available on pip::

    $ pip install ggrecomb

Or you can install ggrecomb from source. See :ref:`Installation<install>` for more information.

In a Python script, import ggrecomb::

    import ggrecomb


2. Make a ParentSequences
-------------------------

Load your parent FASTA file and find additional parents if needed. For this example, we'll use beta-glucosidase (bgl3, PDB ID 1GNX)::

    parent_fn = 'bgl3.fasta'
    p_aln = ggrecomb.ParentSequences.from_fasta(parent_fn)
    p_aln.obtain_seqs(num_final_seqs=4, desired_indentity=0.7)

After running, p_aln is a ParentSequences with four parents that have about 70% pairwise idenity.

See :class:`ggrecomb.ParentSequences` for more options.


3. Run the SCHEMA-RASPP algorithm
---------------------------------

SCHEMA-RASPP finds potential libraries and calculates the probability of Golden Gate assembly for each::

    raspp = ggrecomb.RASPP(p_aln, 5)
    libraries = raspp.vary_m_proxy(60, 100)

This finds libraries with six blocks (five breakpoints) with block sizes between 60 and 100 amino acids.

See ggrecomb.RASPP (TODO: make this link after refactoring PA docstring) for more options.


4. Select and save a library
----------------------------

Let RASPP automatically select a library::

    best_lib = ggrecomb.LibSelector.auto(libraries)
    best_lib.save('bgl3_library_dna.fasta')

The DNA fragments in FASTA format in a file named "bgl3_library_dna.fasta". These fragments are ready to order and assemble with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. You can simulate the Golden Gate reaction using SnapGene.

See ggrecomb.LibSelector for more options.
