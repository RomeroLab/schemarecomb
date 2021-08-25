With ggrecomb, you can design easy-to-use recombinant protein libraries, even if you aren't a computational expert.S

Here's a simple example::

  import ggrecomb

  # Specify library parameters.
  parent_fn = 'P450_sequences.fasta'  # parent_file: 3 parents
  num_blocks = 8  # number of blocks in chimeras
  min_block_len = 40  # min length of chimeric blocks
  max_block_len = 80  # max length of chimeric blocks

  # Create a parent alignment and get the closest PDB structure.
  parents = ggrecomb.ParentSequences.from_fasta(parent_fn, auto_align=True)
  parents.get_PDB()

  # Run SCHEMA-RASPP to get libraries.
  raspp = ggrecomb.RASPP(parents, num_blocks)
  libraries = raspp.vary_m_proxy(min_block_len, max_block_len)

  # Auto-select the best library and save the resulting DNA fragments.
  best_lib = ggrecomb.LibSelector.auto(libraries)
  best_lib.save('library_dna_fragments.fasta')

With this simple script, we generated a three parent, eight block chimeric P450 library. The saved DNA fragments can be ordered directly from a DNA synthesis provider and assembled with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. There's no worrying about adding restriction sites since ggrecomb automatically adds BsaI sites.


Installation
------------

::
    pip install ggrecomb


Documentation
-------------

Package reference material and helpful guides can be found at::
    <insert read the docs link>


Citing
------

If you use ggrecomb in a scientific publication, please cite it as::
    <insert citation format>
