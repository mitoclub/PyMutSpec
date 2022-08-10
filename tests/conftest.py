import pytest

from mutspec.annotation import CodonAnnotation
from mutspec.io import GenomeStates

path_to_states = "./tests/data/states_sample.tsv"
path_to_db = './tests/data/states_sample.db'


@pytest.fixture
def states():
    gs = GenomeStates(path_to_states, path_to_db=path_to_db, mode="dict")
    return gs


@pytest.fixture
def coda():
    gencode = 2
    ca = CodonAnnotation(gencode)
    return ca
