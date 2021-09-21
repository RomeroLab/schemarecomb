With ggrecomb, you can design easy-to-use recombinant protein libraries, even if you aren't a computational expert.

Here's a simple example::

  >>> import ggrecomb as sr
  >>> from Bio import SeqIO
  >>>
  >>> # pytest stuff, you can ignore this.
  >>> getfixture('bgl3_mock_namespace')
  >>> tempdir = getfixture('tmpdir')  #doctest: +ELLIPSIS
  >>> out_fn = tempdir / 'bgl3_dna_frags.fasta'
  >>>
  >>> # Create a parent alignment and get the closest PDB structure.
  >>> fn = 'tests/fixtures/bgl3_1-parent/bgl3_p0.fasta'
  >>> parents = sr.ParentSequences.from_fasta(fn)
  >>> parents.obtain_seqs(6, 0.7)  # BLAST takes about 10 minutes.
  >>> parents.align()  # MUSCLE takes about a minute.
  >>> parents.get_PDB()  # BLAST takes about 10 minutes.
  >>>
  >>> # Run SCHEMA-RASPP to get libraries.
  >>> libraries = sr.generate_libraries(parents, 7)
  >>>
  >>> # Auto-select the best library and save the resulting DNA fragments.
  >>> best_lib = max(libraries, key=lambda x: x.mutation_rate - x.energy)
  >>>
  >>> # Save the generated DNA fragments.
  >>> # out_fn = tempdir + '/' + 'bgl3_dna_frags.fasta'
  >>> SeqIO.write(best_lib.dna_blocks, out_fn, 'fasta')
  42

With this simple script, we generated a six parent, seven block chimeric beta-glucosidase library. The saved DNA fragments can be ordered directly from a DNA synthesis provider and assembled with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. There's no worrying about adding restriction sites since ggrecomb automatically adds BsaI sites.


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
