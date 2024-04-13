import pytest
from ete3 import PhyloTree

from pymutspec.annotation import CodonAnnotation
from pymutspec.io import GenesStates

path_to_states = ["./tests/data/states_sample.tsv"]
path_to_db = "./tests/data/states_sample.db"
path_to_phylip = "./tests/data/example.phy"
path_to_tree_rooted = "./tests/data/treefile_rooted.nwk"


@pytest.fixture
def states():
    gs = GenesStates(path_to_states, mode="dict", use_proba=True)
    return gs


@pytest.fixture
def states_most_probable():
    gs = GenesStates(path_to_states, mode="dict", use_proba=False)
    return gs


@pytest.fixture
def coda():
    gencode = 2
    coda = CodonAnnotation(gencode)
    return coda


@pytest.fixture
def tree_rooted() -> PhyloTree:
    t = PhyloTree(path_to_tree_rooted, format=1)
    return t
