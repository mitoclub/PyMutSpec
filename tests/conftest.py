import pytest
import numpy as np
import pandas as pd

from mutspec.utils.annot import CodonAnnotation
from mutspec.utils.io import GenomeStates

path_to_states = "./data/states_sample.tsv"
path_to_db = './data/states_sample.db'


@pytest.fixture
def states():
    gs = GenomeStates(path_to_states, path_to_db=path_to_db, mode="dict")
    return gs


@pytest.fixture
def coda():
    gencode = 2
    ca = CodonAnnotation(gencode)
    return ca
