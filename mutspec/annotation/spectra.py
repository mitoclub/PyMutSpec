from sys import stderr
from typing import Set, Union, Dict, Iterable

import numpy as np
import pandas as pd

from .mut import CodonAnnotation
from ..constants import possible_sbs12_set, possible_sbs192_set


def calculate_mutspec(
    obs_muts: pd.DataFrame,
    exp_muts_or_genome: Union[Dict[str, float], Iterable],
    gencode: int = None,
    use_context: bool = False,
    use_proba: bool = False,
    verbose=False,
):
    """
    Calculate mutational spectra for mutations dataframe and states frequencies of reference genome

    Arguments
    ---------
    obs_muts: pd.DataFrame
        table containing mutations with annotation; table must contain 2 columns:
        - Mut: str; Pattern: '[ACGT]\[[ACGT]>[ACGT]\][ACGT]'
        - Label: int; [-3, 2]. See CodonAnnotation.get_mut_type
        - ProbaFull (optional, only for use_proba=True) - probability of mutation

    exp_muts_or_genome: dict[str, float] or Iterable
        dictionary that contains expected mutations frequencies of reference genome if use_context=False, 
        else trinucleotide freqs; OR you can pass just genome
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
    _cols = ["Label", "Mut", "ProbaFull"] if use_proba else ["Label", "Mut"]
    for c in _cols:
        assert c in obs_muts.columns, f"Column {c} is not in mut df"

    if isinstance(exp_muts_or_genome, dict):
        exp_muts = exp_muts_or_genome.copy()
    elif isinstance(exp_muts_or_genome, Iterable):
        genome = exp_muts_or_genome.copy()
        if gencode is None:
            raise RuntimeError("If genome passed, gencode argument is required")
        coda = CodonAnnotation(gencode)
        _exp_muts12, _exp_muts192 = coda.collect_exp_mut_freqs(genome)
        exp_muts = _exp_muts192 if use_context else _exp_muts12
    else:
        raise ValueError(
            "'exp_muts_or_genome' must be iterable in case of genome or dict in case of precalculated exp_muts freqs")

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
    mutspec.columns = ["Mut", "ObsFr"]

    # fill unobserved mutations by zeros
    mutspec_appendix = []
    unobserved_sbs = full_sbs.difference(mutspec["Mut"].values)
    for usbs in unobserved_sbs:
        mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
    mutspec = pd.concat([mutspec, pd.DataFrame(mutspec_appendix)], ignore_index=True)

    if use_context:
        sbs = mutspec["Mut"]
        mutspec["Context"] = sbs.str.get(0) + sbs.str.get(2) + sbs.str.get(-1)
    else:
        mutspec["Context"] = mutspec["Mut"].str.get(0)

    mutspec["ExpFr"] = mutspec["Mut"].map(exp_muts)
    mutspec["RawMutSpec"] = (mutspec["ObsFr"] / mutspec["ExpFr"]).fillna(0)
    if verbose:
        for sbs, cnt in mutspec[mutspec.RawMutSpec == np.inf][["Mut", "ObsFr"]].values:
            print(f"WARNING! Substitution {sbs} is unexpected but observed, n={cnt}", file=stderr)
    mutspec["RawMutSpec"] = np.where(mutspec.RawMutSpec == np.inf, mutspec.ObsFr, mutspec.RawMutSpec)
    mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
    mutspec.drop("Context", axis=1, inplace=True)

    assert np.isclose(mutspec.ObsFr.sum(), mut.ProbaFull.sum())
    return mutspec
