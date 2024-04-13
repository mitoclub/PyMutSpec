import numpy as np
import pytest

from pymutspec.io import GenesStates

path_to_states = ["./tests/data/states_sample.tsv"]
path_to_db = './tests/data/states_sample.db'


@pytest.mark.parametrize("use_proba", [True, False])
def test_get_genome_simple(use_proba):
    gs1 = GenesStates(
        path_to_states, path_to_db=path_to_db, 
        mode="db", rewrite=False, use_proba=use_proba,
    )
    gs2 = GenesStates(
        path_to_states, path_to_db=path_to_db, 
        mode="dict", use_proba=use_proba,
    )
    for node in gs1.nodes:
        g1 = gs1.get_genome(node)
        g2 = gs2.get_genome(node)
        for gene in g1:
            assert np.all(g1[gene] == g2[gene])


