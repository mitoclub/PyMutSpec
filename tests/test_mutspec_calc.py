import pytest

import pandas as pd

from mutspec.utils.annot import calculate_mutspec


@pytest.fixture
def mut():
    mut = pd.DataFrame()
    return mut


@pytest.fixture
def freqs():
    fr = dict()
    return 


@pytest.mark.parametrize(
    "lbl",
    [
        pytest.param("all", id='lbl all'),
        pytest.param("syn", id='lbl syn'),
        pytest.param("ff" , id='lbl ff'),
    ]
)
@pytest.mark.parametrize(
    "use_cxt",
    [
        pytest.param(True , id='With cxt'),
        pytest.param(False, id='Without cxt'),
    ]
)
def test_ms_calc_simple(mut, freqs, lbl, use_cxt):
    ms = calculate_mutspec(mut, freqs, lbl, use_cxt, use_proba=False)
    if lbl == ...:
        ...


def test_ms_calc_proba():
    # , use_proba=True
    pass