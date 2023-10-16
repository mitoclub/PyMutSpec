from sys import stderr
from typing import Set, Union, Dict, Iterable

import numpy as np
import pandas as pd

from pymutspec.constants import possible_sbs12_set, possible_sbs192_set, possible_sbs192, possible_sbs12
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
        # TODO add support of sbs12 in Mut column
        mut["Sbs12"] = mut["Mut"].str.slice(2, 5)
        col_mut = "Sbs12"
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

    mutspec["ExpNum"] = mutspec["Mut"].map(exp_muts)
    mutspec["MutSpec"] = (mutspec["ObsNum"] / mutspec["ExpNum"]).fillna(0)
    if verbose:
        msg = mutspec[mutspec["MutSpec"] == np.inf]
        print(f"WARNING! Following substitutions are unexpected but observed:\n{msg}", file=stderr)
    mutspec["MutSpec"] = np.where(mutspec["MutSpec"] == np.inf, 0, mutspec.MutSpec)
    mutspec["MutSpec"] = mutspec["MutSpec"] / mutspec["MutSpec"].sum()
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


def complete_sbs192_columns(df: pd.DataFrame):
    df = df.copy()
    if len(df.columns) != 192:
        for sbs192 in set(possible_sbs192).difference(df.columns.values):
            df[sbs192] = 0.
    df = df[possible_sbs192]
    return df


def jackknife_spectra_sampling(obs: pd.DataFrame, exp: pd.DataFrame, frac=0.5, n=1000):
    if len(obs.columns) == 192 and (obs.columns == possible_sbs192).all() and (exp.columns == possible_sbs192).all():
        assert obs.index.names == ["RefNode", "AltNode"]
        assert exp.index.names == ["Node"]
        altnodes  = obs.index.get_level_values(1).values
        obs_edges = obs
        freqs_nodes = exp
        obs_edges.index = obs_edges.index.reorder_levels(order=["AltNode", "RefNode"])
        freqs_nodes.index.name = "RefNode"
    else:
        altnodes = obs.AltNode.unique()
        obs_edges = obs.groupby(["AltNode", "RefNode", "Mut"]).ProbaFull.sum().unstack()
        obs_edges = complete_sbs192_columns(obs_edges)
        freqs_nodes = exp.rename(columns={"Node": "RefNode"}).groupby(["RefNode", "Mut"]).Proba.sum().unstack()
        freqs_nodes = complete_sbs192_columns(freqs_nodes)

    edges_sample_size = int(len(altnodes) * frac)
    spectra = []
    for _ in range(n):
        altnodes_sample = np.random.choice(altnodes, edges_sample_size, False)
        obs_sample = obs_edges.loc[altnodes_sample].reset_index(0, drop=True)
        exp_sample = freqs_nodes.loc[obs_sample.index]
        
        obs_sample_cnt = obs_sample.sum()
        exp_sample_cnt = exp_sample.sum()

        assert (obs_sample_cnt.index == exp_sample_cnt.index).all()

        sample_spectra = obs_sample_cnt / exp_sample_cnt
        spectra.append(sample_spectra)

    return pd.DataFrame(spectra).fillna(0.)


def collapse_sbs192(df: pd.DataFrame, to=12):
    assert (df.columns == possible_sbs192).all()
    df = df.copy()
    if to == 12:
        for sbs192 in possible_sbs192:
            sbs12 = sbs192[2:5]
            if sbs12 in df.columns.values:
                df[sbs12] += df[sbs192]
            else:
                df[sbs12] = df[sbs192]

        return df[possible_sbs12]
    else:
        raise NotImplementedError()


def calc_edgewise_spectra(obs: pd.DataFrame, exp: pd.DataFrame, nmtypes_cutoff=16, nobs_cuttof=20, collapse_to_12=False, scale=True, both_12_and_192=False):
    if len(obs.columns) == 192 and (obs.columns == possible_sbs192).all() and (exp.columns == possible_sbs192).all():
        assert obs.index.names == ["RefNode", "AltNode"]
        assert exp.index.names == ["Node"]
        obs_edges = obs
        freqs_nodes = exp
        freqs_nodes.index.name = "RefNode"
    else:
        obs_edges = obs.groupby(["RefNode", "AltNode", "Mut"]).ProbaFull.sum().unstack()
        obs_edges = complete_sbs192_columns(obs_edges)
        freqs_nodes = exp.groupby(["Node", "Mut"]).Proba.sum().unstack()
        freqs_nodes.index.name = "RefNode"
        freqs_nodes = complete_sbs192_columns(freqs_nodes)

    if not collapse_to_12:
        obs_edges = obs_edges[((obs_edges > 0).sum(axis=1) >= nmtypes_cutoff) & (obs_edges.sum(axis=1) > nobs_cuttof)]
    
    edges_df = obs_edges.index.to_frame(False)

    freqs_edges = edges_df.merge(freqs_nodes, on="RefNode")\
        .set_index(["RefNode", "AltNode"])[possible_sbs192]
    
    
    assert (obs_edges.columns == freqs_edges.columns).all()
    assert (obs_edges.index == freqs_edges.index).all()

    if both_12_and_192:
        obs_edges12   = collapse_sbs192(obs_edges.fillna(0.),   to=12)
        freqs_edges12 = collapse_sbs192(freqs_edges.fillna(0.), to=12)

        spectra12 = (obs_edges12 / freqs_edges12).replace(np.inf, 0.).fillna(0.)
        spectra192 = (obs_edges / freqs_edges).replace(np.inf, 0.).fillna(0.)
        if scale:
            spectra12 = (spectra12.T / spectra12.T.sum(axis=0)).T
            spectra192 = (spectra192.T / spectra192.T.sum(axis=0)).T

        spectra12 = spectra12.fillna(0)
        spectra192 = spectra192.fillna(0)
        assert not (spectra12 == np.inf).any().any()
        assert not (spectra12.isna()).any().any()
        assert not (spectra192 == np.inf).any().any()
        assert not (spectra192.isna()).any().any()
        
        return spectra12, spectra192

    if collapse_to_12:
        obs_edges   = collapse_sbs192(obs_edges.fillna(0.),   to=12)
        freqs_edges = collapse_sbs192(freqs_edges.fillna(0.), to=12)

    spectra = (obs_edges / freqs_edges).replace(np.inf, 0.).fillna(0.)
    if scale:
        spectra = (spectra.T / spectra.T.sum(axis=0)).T

    spectra = spectra.fillna(0)
    assert not (spectra == np.inf).any().any()
    assert not (spectra.isna()).any().any()
    return spectra


def assign_cat(p: float, interval=0.1):
    if interval < 0.01:
        raise NotImplementedError
    
    left = p // interval / (1 / interval)
    right = left + interval
    if interval >= 0.1:
        return f"{left:.1f}_{right:.1f}"
    else:
        return f"{left:.2f}_{right:.2f}"


def get_cossim(a: pd.DataFrame, b: pd.DataFrame):
    assert (a.columns == b.columns).all()
    
    common_index = a.index.intersection(b.index)
    if len(common_index) == 0:
        return pd.Series()
    
    a = a.loc[common_index]
    b = b.loc[common_index]

    dotprod = (a * b).sum(axis=1)
    a_norm = (a ** 2).sum(axis=1) ** 0.5
    b_norm = (b ** 2).sum(axis=1) ** 0.5
    cossim = dotprod / (a_norm * b_norm)
    return cossim


def get_eucdist(a: pd.DataFrame, b: pd.DataFrame):
    assert (a.columns == b.columns).all()
    
    common_index = a.index.intersection(b.index)
    if len(common_index) == 0:
        return pd.Series()
    
    a = a.loc[common_index]
    b = b.loc[common_index]

    d = ((a - b) ** 2).sum(axis=1) ** 0.5
    return d
