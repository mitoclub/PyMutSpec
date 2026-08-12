import random
import pytest

import pandas as pd

from pymutspec.annotation import calculate_mutspec, calculate_mutrate
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
        divisor = nucl_freqs[lbl].get(sbs, 0)
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
    
    for sbs in cur_mut['Mut'].unique():
        divisor = cxt_freqs[lbl].get(sbs, 0)
        if divisor == 0:
            continue
        cond = cur_mut.Mut.str.fullmatch(sbs.replace("[", r"\[").replace("]", r"\]"))
        if use_proba:
            expected = cur_mut[cond].ProbaFull.sum() / divisor
        else:
            expected = cur_mut[cond].shape[0] / divisor
        observed = ms[ms.Mut == sbs].RawMutSpec.values[0]
        assert observed == expected


def test_calculate_mutspec_keeps_raw_and_scaled(mut, nucl_freqs):
    ms = calculate_mutspec(
        mut, nucl_freqs["all"], use_context=False, use_proba=True,
        drop_underrepresented=False,
    )
    assert "RawMutSpec" in ms.columns
    assert "MutSpec" in ms.columns
    assert pytest.approx(ms["MutSpec"].sum()) == 1.0
    positive = ms[ms["RawMutSpec"] > 0]
    ratios = positive["MutSpec"] / positive["RawMutSpec"]
    assert ratios.max() == pytest.approx(ratios.min())


def test_calculate_mutspec_empty_does_not_nan():
    obs = pd.DataFrame({"Mut": pd.Series(dtype=str), "ProbaFull": pd.Series(dtype=float)})
    exp = {"A>C": 1.0, "C>T": 2.0}
    ms = calculate_mutspec(obs, exp, use_context=False, use_proba=True)
    assert ms["MutSpec"].isna().sum() == 0
    assert (ms["MutSpec"] == 0).all()


def test_calculate_mutspec_accepts_sbs12_mut_column():
    obs = pd.DataFrame({"Mut": ["C>T", "C>T", "A>G"], "ProbaFull": [1.0, 1.0, 0.5]})
    exp = {"C>T": 2.0, "A>G": 1.0}
    ms = calculate_mutspec(obs, exp, use_context=False, use_proba=True, scale=False)
    by_mut = ms.set_index("Mut")
    assert by_mut.loc["C>T", "RawMutSpec"] == pytest.approx(1.0)
    assert by_mut.loc["A>G", "RawMutSpec"] == pytest.approx(0.5)


def test_calculate_mutrate(mut, nucl_freqs):
    rates = calculate_mutrate(mut, nucl_freqs["all"], use_context=False, use_proba=True)
    assert "MutRate" in rates.columns
    assert (rates["MutRate"] == rates["RawMutSpec"]).all()
