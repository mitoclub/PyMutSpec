import numpy as np
import pandas as pd
import pytest

from mutspec import utils
from mutspec.utils import GenomeStates

path_to_states = "./tests/data/states_sample.tsv"
path_to_db = './tests/data/states_sample.db'


@pytest.mark.parametrize(
    "proba_mode",
    [
        pytest.param(False, id='Without proba'),
        pytest.param(True,  id='With proba'),
    ]
)
def test_get_genome_simple(proba_mode):
    gs1 = GenomeStates(
        path_to_states, path_to_db=path_to_db, mode="db", 
        rewrite=False, proba_mode=proba_mode,
    )
    gs2 = GenomeStates(
        path_to_states, path_to_db=path_to_db, mode="dict", 
        proba_mode=proba_mode,
    )
    for node in gs1.nodes:
        g1 = gs1.get_genome(node)
        g2 = gs2.get_genome(node)
        for gene in g1:
            assert np.all(g1[gene] == g2[gene])


