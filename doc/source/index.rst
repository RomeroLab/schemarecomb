.. ggrecomb documentation master file, created by
   sphinx-quickstart on Thu Jul 29 16:15:20 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ggrecomb's documentation!
====================================

With ggrecomb, you can design easy-to-use recombinant protein libraries, even if you aren't a computational expert.


Usage
-----

ggrecomb provides the flexibility to customize recombinant library design with a number of options. Check out the :ref:`Biologist's Guide<biologists>` if you're comfortable with recombinant protein library design but want to be guided through using the library. Conversely, users with computational expertise but limited biological experience should read the :ref:`Programmer's Guide<programmers>`. Finally, if you know what you're doing, glance at the :ref:`Quickstart <quickstart>` document and look at the :ref:`Reference Manual<reference>` as needed.

Here's a simple example::

  import ggrecomb

  # Specify library parameters.
  parent_fn = 'P450_sequences.fasta'  # parent_file: 3 parents
  num_blocks = 8  # number of blocks in chimeras
  min_block_len = 40  # min length of chimeric blocks
  max_block_len = 80  # max length of chimeric blocks

  # Create a parent alignment.
  p_aln = ggrecomb.ParentAlignment.from_fasta(parent_fn)

  # Run SCHEMA-RASPP to get libraries.
  raspp = ggrecomb.RASPP(p_aln, num_blocks)
  libraries = raspp.vary_m_proxy(min_block_len, max_block_len)

  # Auto-select the best library and save the resulting DNA fragments.
  best_lib = ggrecomb.LibSelector.auto(libraries)
  best_lib.save('library_dna_fragments.fasta')

With this simple script, we generated a three parent, eight block chimeric P450 library. The saved DNA fragments can be ordered directly from a DNA synthesis provider and assembled with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. There's no worrying about adding restriction sites since ggrecomb automatically adds BsaI sites.


Contents
--------

.. toctree::
   :maxdepth: 2
    
   quickstart
   install
   biologists
   programmers
   reference


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
