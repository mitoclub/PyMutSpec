from sys import stderr
from typing import Set, Union, Dict, Iterable

import numpy as np
import pandas as pd

from pymutspec.constants import possible_sbs12_set, possible_sbs192_set
from .mut import CodonAnnotation
from .auxiliary import rev_comp


def calculate_mutspec(
    obs_muts: pd.DataFrame,
    exp_muts: Dict[str, float],
    use_context: bool = False,
    use_proba: bool = False,
    verbose=False,
    fill_unobserved=True,
):
    """
    Calculate mutational spectra for mutations dataframe and states frequencies of reference genome

    Arguments
    ---------
    obs_muts: pd.DataFrame
        table containing mutations with annotation; table must contain 2 columns:
        - Mut: str; Pattern: '[ACGT]\[[ACGT]>[ACGT]\][ACGT]'
        - ProbaFull (optional, only for use_proba=True) - probability of mutation

    exp_muts: dict[str, float]
        dictionary that contains expected mutations frequencies of reference genome if use_context=False, 
        else trinucleotide freqs
    label: str
        kind of needed mutspec, coulb be one of ['all', 'syn', 'ff']
    gencode: int
        Number of genetic code to use in expected mutations collection, required if exp_muts_or_genome is genome
    use_context: bool
        To use trinucleotide context or not, in other words calculate 192 component mutspec
    use_proba: bool
        To use probabilities of mutations or not. Usefull if you have such probabiliies
    verbose: bool
        Show warning messages or not

    Return
    -------
    mutspec: pd.DataFrame
        table, containing extended mutspec values including observed mutations numbers. 
        If use_context=True len(mutspec) = 192, else len(mutspec) = 12
    """
    _cols = ["Mut", "ProbaFull"] if use_proba else ["Mut"]
    for c in _cols:
        assert c in obs_muts.columns, f"Column {c} is not in mut df"

    if not isinstance(exp_muts, dict):
        raise ValueError("'exp_muts' must be dict with mutations freqs")

    mut = obs_muts.copy()
    if use_context:
        col_mut = "Mut"
        full_sbs = possible_sbs192_set
    else:
        mut["MutBase"] = mut["Mut"].str.slice(2, 5)
        col_mut = "MutBase"
        full_sbs = possible_sbs12_set

    if not use_proba:
        mut["ProbaFull"] = 1

    mutspec = mut.groupby(col_mut)["ProbaFull"].sum().reset_index()
    mutspec.columns = ["Mut", "ObsNum"]

    if fill_unobserved:
        # fill unobserved mutations by zeros
        mutspec_appendix = []
        unobserved_sbs = full_sbs.difference(mutspec["Mut"].values)
        for usbs in unobserved_sbs:
            mutspec_appendix.append({"Mut": usbs, "ObsNum": 0})
        mutspec = pd.concat([mutspec, pd.DataFrame(mutspec_appendix)], ignore_index=True)

    if use_context:
        sbs = mutspec["Mut"]
        mutspec["Context"] = sbs.str.get(0) + sbs.str.get(2) + sbs.str.get(-1)
    else:
        mutspec["Context"] = mutspec["Mut"].str.get(0)

    mutspec["ExpNum"] = mutspec["Mut"].map(exp_muts)
    mutspec["RawMutSpec"] = (mutspec["ObsNum"] / mutspec["ExpNum"]).fillna(0)
    if verbose:
        for sbs, cnt in mutspec[mutspec.RawMutSpec == np.inf][["Mut", "ObsNum"]].values:
            print(f"WARNING! Substitution {sbs} is unexpected but observed, n={cnt}", file=stderr)
    mutspec["RawMutSpec"] = np.where(mutspec.RawMutSpec == np.inf, mutspec.ObsNum, mutspec.RawMutSpec) # TODO test 
    mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
    mutspec.drop("Context", axis=1, inplace=True)

    assert np.isclose(mutspec.ObsNum.sum(), mut.ProbaFull.sum())
    return mutspec


def collapse_mutspec(ms192: pd.DataFrame):
    assert ms192.shape[0] == 192
    for c in ["Mut", "ObsFr", "ExpFr"]:
        assert c in ms192.columns

    ms1 = ms192[ms192["Mut"].str.get(2).isin(list("CT"))]
    ms2 = ms192[ms192["Mut"].str.get(2).isin(list("AG"))]
    ms2["Mut"] = ms2["Mut"].apply(rev_comp)

    ms96 = pd.concat([ms1, ms2]).groupby("Mut")[["ObsFr", "ExpFr"]].sum()
    ms96["RawMutSpec"] = ms96["ObsFr"] / ms96["ExpFr"]
    ms96["MutSpec"] = ms96["RawMutSpec"] / ms96["RawMutSpec"].sum()
    ms96 = ms96.fillna(0).replace(np.inf, 0)
    return ms96
