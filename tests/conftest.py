import os
from pathlib import Path

import pytest

fixture_dir = Path(os.path.dirname(__file__)) / 'fixtures/'


@pytest.fixture()
def bgl3_fasta_filename():
    return fixture_dir / 'bgl3_full' / 'bgl3_sequences.fasta'


@pytest.fixture()
def bgl3_pdb_filename():
    return fixture_dir / 'bgl3_full' / '1GNX.pdb'
