.. schemarecomb documentation master file, created by
   sphinx-quickstart on Thu Jul 29 16:15:20 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to schemarecomb's documentation!
========================================

With schemarecomb, you can design easy-to-use recombinant protein libraries, even if you aren't a computational expert.


Usage
-----

..
  schemarecomb provides the flexibility to customize recombinant library design with a number of options. Check out the :ref:`Biologist's Guide<biologists>` if you're comfortable with recombinant protein library design but want to be guided through using the library. Conversely, users with computational expertise but limited biological experience should read the :ref:`Programmer's Guide<programmers>`. Finally, if you know what you're doing, glance at the :ref:`Quickstart <quickstart>` document and look at the :ref:`Reference Manual<reference>` as needed.


Here's a simple example::

  >>> import schemarecomb as sr
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
  >>> SeqIO.write(best_lib.dna_blocks, out_fn, 'fasta')
  42

With this simple script, we generated a six parent, seven block chimeric beta-glucosidase library. The saved DNA fragments can be ordered directly from a DNA synthesis provider and assembled with `NEB's Golden Gate Assembly Kit <https://www.neb.com/products/e1601-neb-golden-gate-assembly-mix>`_. There's no worrying about adding restriction sites since schemarecomb automatically adds BsaI sites.

View the :ref:`Quickstart Guide <quickstart>` for more example scripts and the :ref:`Reference Manual <reference>` for more details on specific classes and modules.

Biologist's, Programmer's, and Installation guides coming in version 0.2.0.


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
