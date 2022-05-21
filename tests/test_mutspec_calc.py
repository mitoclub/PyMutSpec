import random
import pytest

import pandas as pd

from mutspec.utils.annot import calculate_mutspec
from mutspec.utils.constants import possible_sbs12, possible_sbs192, possible_codons


@pytest.fixture
def mut():
    data = [
        [0, "A[A>T]G", 0.13],
        [1, "C[C>G]T", 0.06],
        [2, "G[C>T]C", 0.62],
        [0, "C[C>G]A", 0.83],
        [1, "T[A>G]C", 0.37],
        [0, "G[G>T]G", 0.01],
        [2, "A[C>T]A", 0.91],
        [1, "T[C>A]C", 0.27],
        [0, "C[T>A]T", 0.45],
        [2, "T[G>T]G", 0.99],
    ]
    mut = pd.DataFrame(data, columns=["Label", "Mut", "ProbaFull"])
    return mut


@pytest.fixture
def nucl_freqs():
    fr = {
        "all": {"A": 2, "C": 8, "G": 4, "T": 3},
        "syn": {"A": 1, "C": 5, "G": 3, "T": 2},
        "ff" : {"C": 3, "G": 3, "T": 1},
    }
    return fr


@pytest.fixture
def cxt_freqs():
    fr = dict()
    for lbl, max_num in zip(["all", "syn", "ff"], [20, 13, 7]):
        fr[lbl] = {cxt: max(0, random.randint(-5, max_num)) for cxt in possible_codons}
    return fr


@pytest.mark.parametrize(
    "use_proba",
    [
        pytest.param(True , id='With proba'),
        pytest.param(False, id='Without proba'),
    ]
)
def test_ms12_calc(mut, nucl_freqs, use_proba):
    for lbl_id, lbl in enumerate(("all", "syn", 'ff')):
        ms = calculate_mutspec(
            mut, nucl_freqs[lbl], lbl, 
            use_context=False, use_proba=use_proba,
        )
        cur_mut = mut[(mut.Label >= lbl_id)]
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


@pytest.mark.parametrize(
    "use_proba",
    [
        pytest.param(True , id='With proba'),
        pytest.param(False, id='Without proba'),
    ]
)
def test_ms192_calc(mut, cxt_freqs, use_proba):
    # TODO
    for lbl_id, lbl in enumerate(("all", "syn", 'ff')):
        ms = calculate_mutspec(
            mut, cxt_freqs[lbl], lbl, 
            use_context=False, use_proba=use_proba,
        )
        cur_mut = mut[(mut.Label >= lbl_id)]
        for sbs in possible_sbs192:
            cxt = sbs[0] + sbs[2] + sbs[-1]
            divisor = cxt_freqs[lbl].get(cxt, 0)
            if divisor <= 0:
                continue
            if use_proba:
                expected = cur_mut[(cur_mut.Mut.str.contains(sbs))].ProbaFull.sum() / divisor
            else:
                expected = cur_mut[(cur_mut.Mut.str.contains(sbs))].shape[0] / divisor
            observed = ms[ms.Mut == sbs].RawMutSpec.values[0]
            assert observed == expected
