import random
import pytest

import pandas as pd

from pymutspec.annotation import calculate_mutspec
from pymutspec.constants import possible_sbs12, possible_sbs192


@pytest.fixture
def mut():
    data = [
        [0, "A[A>T]G", 0.93],
        [1, "C[C>G]T", 0.86],
        [2, "G[C>T]C", 0.42],
        [0, "C[C>G]A", 0.63],
        [1, "T[A>G]C", 0.97],
        [0, "G[G>T]G", 0.41],
        [2, "A[C>T]A", 0.91],
        [1, "T[C>A]C", 0.57],
        [0, "C[T>A]T", 0.65],
        [2, "T[G>T]G", 0.39],
    ]
    mut = pd.DataFrame(data, columns=["Label", "Mut", "ProbaFull"])
    return mut


@pytest.fixture
def nucl_freqs():
    fr = {
        "all": {"A>C": 2, "A>T": 1, "C>T": 8, "C>A": 2, "C>G": 1, "G>A": 4, "G>T": 2, "T>C": 3, "T>A": 1},
        "syn": {"A>C": 1, "C>T": 6, "C>A": 1, "G>A": 1, "G>T": 1, "T>C": 2},
        "ff" : {"A>C": 1, "C>T": 3, "C>A": 1, "G>T": 1},
    }
    return fr


@pytest.fixture
def cxt_freqs():
    fr = dict()
    for lbl, max_num in zip(["all", "syn", "ff"], [20, 13, 7]):
        fr[lbl] = {cxt: max(0, random.randint(-5, max_num)) for cxt in possible_sbs192}
    return fr


@pytest.mark.parametrize("use_proba", [False, True])
@pytest.mark.parametrize("lbl_id", [0, 1, 2])
def test_ms12_calc(mut, nucl_freqs, use_proba, lbl_id):
    """test only RawMutSpec values"""
    if lbl_id == 0:
        lbl = "all"
    elif lbl_id == 1:
        lbl = "syn"
    elif lbl_id == 2:
        lbl = "ff"
    cur_mut = mut[(mut.Label >= lbl_id)]
    ms = calculate_mutspec(cur_mut, nucl_freqs[lbl], use_context=False, use_proba=use_proba)
    for sbs in possible_sbs12:
        divisor = nucl_freqs[lbl].get(sbs[0], 0)
        if divisor <= 0:
            continue
        if use_proba:
            expected = cur_mut[(cur_mut.Mut.str.contains(sbs))].ProbaFull.sum() / divisor
        else:
            expected = cur_mut[(cur_mut.Mut.str.contains(sbs))].shape[0] / divisor
        observed = ms[ms.Mut == sbs].RawMutSpec.values[0]
        assert observed == expected        


@pytest.mark.parametrize("use_proba", [True, False])
@pytest.mark.parametrize("lbl_id", [0, 1, 2])
def test_ms192_calc(mut, cxt_freqs, use_proba, lbl_id):
    """test only RawMutSpec values"""
    if lbl_id == 0:
        lbl = "all"
    elif lbl_id == 1:
        lbl = "syn"
    elif lbl_id == 2:
        lbl = "ff"
    cur_mut = mut[(mut.Label >= lbl_id)]
    ms = calculate_mutspec(cur_mut, cxt_freqs[lbl], use_context=True, 
                           use_proba=use_proba, fill_unobserved=False)
    
    for sbs in mut['Mut'].unique():
        cxt = sbs[0] + sbs[2] + sbs[-1]
        divisor = cxt_freqs[lbl].get(cxt, 0)
        if divisor == 0:
            continue
        cond = cur_mut.Mut.str.fullmatch(sbs.replace("[", "\[").replace("]", "\]"))
        if use_proba:
            expected = cur_mut[cond].ProbaFull.sum() / divisor
        else:
            expected = cur_mut[cond].shape[0] / divisor
        observed = ms[ms.Mut == sbs].RawMutSpec.values[0]
        assert observed == expected
