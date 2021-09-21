import schemarecomb as sr


def test_defaults(bgl3_parent_alignment):
    libraries = sr.generate_libraries(bgl3_parent_alignment, 8)

    print()
    print(len(libraries))


"""
def test_raspp(bgl3_mock_namespace):
    # Specify library parameters.
    from Bio import SeqIO
    fn = 'tests/fixtures/bgl3_1-parent/bgl3_p0.fasta'

    # Create a parent alignment and get the closest PDB structure.
    print('before setup')
    parents = sr.ParentSequences.from_fasta(fn)
    parents.obtain_seqs(6, 0.7)  # BLAST takes about 10 minutes.
    parents.align()  # MUSCLE takes about a minute.
    parents.get_PDB()  # BLAST takes about 10 minutes.

    print('before generate')

    # Run SCHEMA-RASPP to get libraries.
    libraries = sr.generate_libraries(parents, 7)

    print('selection')

    # Auto-select the best library and save the resulting DNA fragments.
    best_lib = max(libraries, key=lambda x: x.mutation_rate - x.energy)
    SeqIO.write(best_lib.dna_blocks, 'bgl3_dna_frags.fasta', 'fasta')
"""
