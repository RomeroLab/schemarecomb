import os
from pathlib import Path

from Bio import SeqIO
import pytest


@pytest.fixture
def fixture_dir():
    return Path(os.path.dirname(__file__)) / 'tests' / 'fixtures/'


@pytest.fixture
def bgl3_fasta_filename(fixture_dir):
    return fixture_dir / 'bgl3_full' / 'bgl3_sequences.fasta'


@pytest.fixture
def bgl3_records(bgl3_fasta_filename):
    return list(SeqIO.parse(bgl3_fasta_filename, 'fasta'))


@pytest.fixture
def bgl3_pdb_filename(fixture_dir):
    return fixture_dir / 'bgl3_full' / '1GNX.pdb'
