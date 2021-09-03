With ggrecomb, you can design easy-to-use recombinant protein libraries, even if you aren't a computational expert.

Here's a simple example::

  >>> import ggrecomb
  >>>
  >>> # Specify library parameters.
  >>> parent_fn = 'P450_sequences.fasta'  # parent_file: 3 parents
  >>> num_blocks = 8  # number of blocks in chimeras
  >>> min_block_len = 40  # min length of chimeric blocks
  >>> max_block_len = 80  # max length of chimeric blocks
  >>>
  >>> # Create a parent alignment and get the closest PDB structure.
  >>> parents = ggrecomb.ParentSequences.from_fasta(parent_fn, auto_align=True)
  >>> parents.get_PDB()
  >>>
  >>> # Run SCHEMA-RASPP to get libraries.
  >>> raspp = ggrecomb.RASPP(parents, num_blocks)
  >>> library_list = raspp.vary_m_proxy(min_block_len, max_block_len)
  >>>
  >>> # Auto-select the best library and save the resulting DNA fragments.
  >>> best_lib = ggrecomb.libraries.auto_select(library_list)
  >>> best_lib.save('library_dna_fragments.fasta')

With this simple script, we generated a three parent, eight block chimeric P450 library. The saved DNA fragments can be ordered directly from a DNA synthesis provider and assembled with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. There's no worrying about adding restriction sites since ggrecomb automatically adds BsaI sites.


Why Recombinant Proteins?
-------------------------

Engineering proteins with recombinant libraries has a number of advantages over traditional directed evolution. Primarily, the use of multiple parental sequences results in a more "global" optimization over the protein function space. That is, recombinant libraries are less limited by local optima and environmental dependence when compared to single-gene mutagenesis. Furthermore, the recombinant libraries generally have a greater fraction of functional variants at the same level of diversity. Recombinant variants may also be assembled in isolation if desired, without requiring a transformation. The relative advantages of single-gene mutagenesis, such as fine-grain sequence space exploration, may be mitigated by subjecting a recombinant library to the same mutagenesis techniques.

So why doesn't everybody do recombinant protein engineering? Historically, there's a number of technical challenges that made recombinant libraries impractical for general use. Namely, a great deal of computational expertise and time is needed to manually generate and select suitable libraries. Even with the required computational resources, mutagenesis was significantly easier than assembling 8+ DNA fragments when protein engineering was developing its fundamentals, so nearly everybody opted for traditional directed evolution and passed that practice down to their students and mentees.

The goal of this software package is to make recombinant library design accessible and convenient for protein engineers of all computational skill levels. ggrecomb designs libraries that come ready to order and construct with a simple Golden Gate Assembly reaction. To learn more, read the `ggrecomb documentation <https://ggrecomb.readthedocs.io/en/latest/>`_.


Installation
------------

.. code-block:: bash

    $ pip install ggrecomb


Documentation
-------------

Package reference material and helpful guides can be found at:

    https://ggrecomb.readthedocs.io/en/latest/


Citing
------

..
    https://www.software.ac.uk/how-cite-software?_ga=1.54830891.1882560887.1489012280

If you use ggrecomb in a scientific publication, please cite it as::

    Bremer, B. & Romero, P. (2021). ggrecomb [Software]. Available from https://github.com/RomeroLab/ggrecomb.
