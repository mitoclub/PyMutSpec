import numpy as np
import pandas as pd

from mutspec import utils
from mutspec.utils import GenomeStates

path_to_states = "./data/states_sample.tsv"
path_to_db = './data/states_sample.db'


def test_get_genome_correct():
    gs1 = GenomeStates(
        path_to_states, path_to_db=path_to_db, mode="db", 
        rewrite=False,
    )
    gs2 = GenomeStates(
        path_to_states, path_to_db=path_to_db, mode="dict"
    )
    for node in gs1.nodes:
        g1 = gs1.get_genome(node)
        g2 = gs2.get_genome(node)
        for gene in g1:
            assert np.all(g1[gene] == g2[gene])
