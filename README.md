# SCHEMA-library-design
Tools for designing site-directed recombination libraries using SCHEMA and RASPP

General Workflow:
  1.  Generate pickled multiple sequence alignment (MSA) and contacts files. The MSA file should be a 2D list where every element of       the outer list contains the amino acids of all sequences at a specific position. The contacts file numbering must be indexed from 0 and correspond to the MSA positional numbering.
  2.  Use run_raspp.py with the pickled MSA and contacts files as input to generate SCHEMA_gg libraries.
  3.  Run read_raspp_GG_output.py to select a library and save it as a pickled file.
  4.  Use generate_gg_dna.py with the codon-aligned nucleotide sequences to generate golden gate fragments for ordering.

_Forked from RomeroLab for personal modification. Will probably merge later._
