from sys import stderr
from typing import Dict

import numpy as np
import pandas as pd

from ..constants import (
    possible_sbs12_set, possible_sbs192_set, 
    possible_sbs192, possible_sbs12
)
from .auxiliary import rev_comp


_SBS12_RE = r"[ACGT]>[ACGT]"


def _as_sbs12(mut_series: pd.Series) -> pd.Series:
    """Accept both 192-component (``A[C>T]G``) and 12-component (``C>T``) Mut values."""
    mut_str = mut_series.astype(str)
    sbs12 = mut_str.str.slice(2, 5)
    already_sbs12 = mut_str.str.fullmatch(_SBS12_RE)
    sbs12 = sbs12.where(~already_sbs12, mut_str)
    return sbs12


def calculate_mutspec(
    obs_muts: pd.DataFrame,
    exp_muts: Dict[str, float],
    use_context: bool = False,
    use_proba: bool = False,
    scale=True,
    fill_unobserved=True,
    drop_underrepresented=True,
    nobs_min=0.9,
    nexp_min=0.9,
    verbose=False,
):
    """
    Calculate mutational spectra for mutations dataframe and states frequencies of reference genome

    Arguments
    ---------
    obs_muts: pd.DataFrame
        table containing mutations with annotation; table must contain 2 columns:
        - Mut: str; Pattern: ``[ACGT]\\[[ACGT]>[ACGT]\\][ACGT]`` or ``[ACGT]>[ACGT]``
        - ProbaFull (optional, only for use_proba=True) - probability of mutation

    exp_muts: dict[str, float]
        dictionary that contains expected mutations frequencies of reference genome if use_context=False, 
        else trinucleotide freqs
    use_context: bool
        To use trinucleotide context or not, in other words calculate 192 component mutspec
    use_proba: bool
        To use probabilities of mutations or not. Usefull if you have such probabiliies
    scale: bool
        Scale spectrum vector (divide by sum)
    fill_unobserved: bool
        Fill table with mutation types that didn't observed
    drop_underrepresented: bool
        Drop underrepresented mutation types from spectrum according to `nobs_min` and `nexp_min`
    nobs_min: float/int
        Minimal number of observed mutations for each mutation type
    nexp_min: float/int
        Minimal number of expected mutations for each mutation type
    verbose: bool
        Show warning messages

    Return
    -------
    mutspec: pd.DataFrame
        table, containing observed/expected counts and both unscaled (``RawMutSpec``)
        and optionally scaled (``MutSpec``) rates.
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
        mut["Sbs12"] = _as_sbs12(mut["Mut"])
        col_mut = "Sbs12"
        full_sbs = possible_sbs12_set

    if not use_proba:
        mut["ProbaFull"] = 1

    mutspec = mut.groupby(col_mut, sort=False)["ProbaFull"].sum().reset_index()
    mutspec.columns = ["Mut", "ObsNum"]

    if fill_unobserved:
        unobserved_sbs = full_sbs.difference(mutspec["Mut"].values)
        if unobserved_sbs:
            mutspec = pd.concat(
                [mutspec, pd.DataFrame({"Mut": list(unobserved_sbs), "ObsNum": 0})],
                ignore_index=True,
            )

    mutspec["ExpNum"] = mutspec["Mut"].map(exp_muts)
    raw = mutspec["ObsNum"] / mutspec["ExpNum"]
    if verbose:
        unexpected = mutspec[(mutspec["ObsNum"] > 0) & (mutspec["ExpNum"].fillna(0) <= 0)]
        if len(unexpected) > 0:
            print(f"WARNING! Following substitutions are unexpected but observed:\n{unexpected}", file=stderr)

    raw = raw.replace([np.inf, -np.inf], np.nan).fillna(0)
    mutspec["RawMutSpec"] = raw
    mutspec["MutSpec"] = raw

    if drop_underrepresented:
        under = (mutspec["ObsNum"] < nobs_min) | (mutspec["ExpNum"] < nexp_min)
        mutspec.loc[under.fillna(False), "MutSpec"] = 0.
    if scale:
        total = mutspec["MutSpec"].sum()
        if total > 0:
            mutspec["MutSpec"] = mutspec["MutSpec"] / total
        else:
            mutspec["MutSpec"] = 0.

    return mutspec


def calculate_mutrate(
    obs_muts: pd.DataFrame,
    exp_muts: Dict[str, float],
    use_context: bool = False,
    use_proba: bool = False,
    fill_unobserved=True,
    verbose=False,
):
    """
    Calculate mutation rates (observed / expected) without scaling to a spectrum.

    This is the un-normalised rate vector that ``calculate_mutspec`` stores in
    ``RawMutSpec``.  The returned table uses the column name ``MutRate``.

    Arguments
    ---------
    obs_muts, exp_muts, use_context, use_proba, fill_unobserved, verbose
        Same meaning as in :func:`calculate_mutspec`.

    Return
    ------
    mutrate: pd.DataFrame
        Table with columns ``Mut``, ``ObsNum``, ``ExpNum``, ``RawMutSpec``,
        ``MutSpec`` (unscaled) and ``MutRate`` (alias of ``RawMutSpec``).
    """
    mutrate = calculate_mutspec(
        obs_muts,
        exp_muts,
        use_context=use_context,
        use_proba=use_proba,
        scale=False,
        fill_unobserved=fill_unobserved,
        drop_underrepresented=False,
        verbose=verbose,
    )
    mutrate["MutRate"] = mutrate["RawMutSpec"]
    return mutrate


def sample_spectrum(obs_df: pd.DataFrame, exp_freqs,
                    use_proba=True, use_context=False, 
                    frac=0.5, nreplics=100):
    """
    Sample half of branches and calculate spectrum
    """
    samples = []
    edges = obs_df.AltNode.unique()
    n_to_sample = int(len(edges) * frac)
    for _ in range(nreplics):
        cur_edges = np.random.choice(edges, n_to_sample, replace=False)
        obs_smpl = obs_df[obs_df['AltNode'].isin(cur_edges)]
        one_spectrum = calculate_mutspec(
            obs_smpl, exp_freqs, use_context=use_context, use_proba=use_proba)
        samples.append(one_spectrum)

    sampled = pd.concat(samples)

    quartiles = sampled.groupby('Mut')['MutSpec'].quantile([0.05, 0.5, 0.95]).unstack().rename(
        columns={0.05: "MutSpec_q05", 0.5: "MutSpec_median", 0.95: "MutSpec_q95"}).reset_index()
    return quartiles
    

def get_iqr_bounds(series: pd.Series):
    "Function calculates Interquartile range (IQR) and used for outliers filtrtation"
    q1 = series.quantile(0.25)
    q3 = series.quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return lower_bound, upper_bound


def filter_outlier_branches(obs_df: pd.DataFrame, use_proba=True):
    """
    Remove branches with an outlier-high number of observed mutations.

    Outliers are identified using the IQR method: branches whose mutation
    count exceeds ``Q3 + 1.5 * IQR`` are dropped.

    Arguments
    ---------
    obs_df: pd.DataFrame
        Observed-mutations table containing at least the columns
        ``'AltNode'``, ``'Mut'``, and optionally ``'ProbaMut'`` or
        ``'ProbaFull'``.
    use_proba: bool
        If ``True`` sum ``'ProbaMut'`` (falling back to ``'ProbaFull'``)
        per branch; otherwise count rows.

    Return
    ------
    obs_df_flt: pd.DataFrame
        Filtered mutations table with outlier branches removed.
    """
    if use_proba:
        if "ProbaMut" in obs_df.columns:
            proba_col = "ProbaMut"
        elif "ProbaFull" in obs_df.columns:
            proba_col = "ProbaFull"
        else:
            raise ValueError(
                "use_proba=True requires a 'ProbaMut' or 'ProbaFull' column"
            )
        edge_nobs = obs_df.groupby('AltNode')[proba_col].sum()
    else:
        edge_nobs = obs_df.groupby('AltNode')['Mut'].count()

    _, upper_bound = get_iqr_bounds(edge_nobs)
    edge_nobs_flt = edge_nobs[edge_nobs < upper_bound]
    selected_branches = edge_nobs_flt.index
    obs_df_flt = obs_df[obs_df['AltNode'].isin(selected_branches)]
    return obs_df_flt


def collapse_mutspec(ms192: pd.DataFrame):
    """
    Collapse a 192-component spectrum to 96 components using reverse complement.

    Mutations on the ``A``/``G`` strand are reverse-complemented so that all
    substitutions are expressed relative to the pyrimidine base (``C`` or
    ``T``), then the ``ObsFr`` and ``ExpFr`` columns are summed for matching
    contexts, yielding a 96-component spectrum.

    Arguments
    ---------
    ms192: pd.DataFrame
        192-component spectrum table.  Must contain column ``'Mut'`` and
        either ``'ObsFr'``/``'ExpFr'`` or ``'ObsNum'``/``'ExpNum'``, and
        must have exactly 192 rows.

    Return
    ------
    ms96: pd.DataFrame
        96-component spectrum with columns ``'ObsFr'``, ``'ExpFr'``,
        ``'RawMutSpec'``, and ``'MutSpec'`` (normalised to sum to 1).

    Raises
    ------
    AssertionError
        If ``ms192`` does not have exactly 192 rows or is missing required
        columns.
    """
    assert ms192.shape[0] == 192, f"Expected 192 rows, got {ms192.shape[0]}"
    ms192 = ms192.copy()
    if "ObsFr" not in ms192.columns and "ObsNum" in ms192.columns:
        ms192["ObsFr"] = ms192["ObsNum"]
    if "ExpFr" not in ms192.columns and "ExpNum" in ms192.columns:
        ms192["ExpFr"] = ms192["ExpNum"]
    for c in ["Mut", "ObsFr", "ExpFr"]:
        assert c in ms192.columns, f"Required column '{c}' not found in ms192"

    ms1 = ms192.loc[ms192["Mut"].str.get(2).isin(list("CT"))].copy()
    ms2 = ms192.loc[ms192["Mut"].str.get(2).isin(list("AG"))].copy()
    ms2["Mut"] = ms2["Mut"].apply(rev_comp)

    ms96 = pd.concat([ms1, ms2]).groupby("Mut")[["ObsFr", "ExpFr"]].sum()
    ms96["RawMutSpec"] = ms96["ObsFr"] / ms96["ExpFr"]
    ms96["MutSpec"] = ms96["RawMutSpec"] / ms96["RawMutSpec"].sum()
    ms96 = ms96.fillna(0).replace(np.inf, 0)
    return ms96


def complete_sbs192_columns(df: pd.DataFrame):
    """
    Ensure a DataFrame has all 192 SBS columns, filling missing ones with 0.

    The resulting DataFrame is reordered so its columns follow the canonical
    ``possible_sbs192`` order.

    Arguments
    ---------
    df: pd.DataFrame
        DataFrame whose columns are a (possibly incomplete) subset of the 192
        SBS mutation types.

    Return
    ------
    df: pd.DataFrame
        DataFrame with exactly 192 columns in canonical order.
    """
    df = df.copy()
    missing = [sbs for sbs in possible_sbs192 if sbs not in df.columns]
    if missing:
        df = pd.concat(
            [df, pd.DataFrame(0.0, index=df.index, columns=missing)],
            axis=1,
        )
    return df[possible_sbs192]


def collapse_sbs192(df: pd.DataFrame, to=12):
    """
    Sum a 192-component SBS DataFrame into a 12-component representation.

    Each 192-component mutation type is mapped to its 12-component base
    substitution (the middle three characters, e.g. ``'C>A'``), and the
    values are accumulated.

    Arguments
    ---------
    df: pd.DataFrame
        DataFrame with columns equal to ``possible_sbs192`` in canonical order.
        Each row typically represents one sample or branch.
    to: int
        Target number of components.  Currently only ``12`` is supported.

    Return
    ------
    df12: pd.DataFrame
        DataFrame with 12 columns corresponding to the 12 base substitution
        types in ``possible_sbs12`` order.

    Raises
    ------
    AssertionError
        If ``df.columns`` does not match ``possible_sbs192``.
    NotImplementedError
        If ``to`` is not ``12``.
    """
    assert (df.columns == possible_sbs192).all(), \
        "DataFrame columns must match possible_sbs192 in canonical order"
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


def jackknife_spectra_sampling(obs: pd.DataFrame, exp: pd.DataFrame, frac=0.5, n=1000):
    """
    Estimate spectrum variability via jackknife resampling of tree branches.

    On each iteration a random subset of branches (edges) is drawn without
    replacement and a per-branch spectrum ratio ``obs / exp`` is computed.
    The resulting collection of spectra can be used to derive confidence
    intervals.

    Arguments
    ---------
    obs: pd.DataFrame
        Observed mutations.  Either a pre-pivoted wide DataFrame with 192
        SBS columns and a ``(RefNode, AltNode)`` MultiIndex, or a long-format
        DataFrame with columns ``'AltNode'``, ``'RefNode'``, ``'Mut'``,
        and ``'ProbaFull'``.
    exp: pd.DataFrame
        Expected mutation frequencies.  Either a pre-pivoted wide DataFrame
        with 192 SBS columns and a ``Node`` index, or a long-format DataFrame
        with columns ``'Node'``, ``'Mut'``, and ``'Proba'``.
    frac: float
        Fraction of branches to sample on each iteration.
    n: int
        Number of jackknife iterations.

    Return
    ------
    spectra: pd.DataFrame
        DataFrame of shape ``(n, 192)`` where each row is the spectrum
        computed from one jackknife sample.
    """
    if len(obs.columns) == 192 and \
            (obs.columns == possible_sbs192).all() and \
                (exp.columns == possible_sbs192).all():
        assert obs.index.names == ["RefNode", "AltNode"]
        assert exp.index.names == ["Node"]
        altnodes  = obs.index.get_level_values(1).values
        obs_edges = obs.copy()
        freqs_nodes = exp.copy()
        obs_edges.index = obs_edges.index.reorder_levels(order=["AltNode", "RefNode"])
        freqs_nodes.index.name = "RefNode"
    else:
        altnodes = obs.AltNode.unique()
        obs_edges = obs.groupby(["AltNode", "RefNode", "Mut"]).ProbaFull.sum().unstack()
        obs_edges = complete_sbs192_columns(obs_edges)
        freqs_nodes = exp.rename(columns={"Node": "RefNode"})\
            .groupby(["RefNode", "Mut"]).Proba.sum().unstack()
        freqs_nodes = complete_sbs192_columns(freqs_nodes)

    edges_sample_size = int(len(altnodes) * frac)
    spectra = []
    for _ in range(n):
        altnodes_sample = np.random.choice(altnodes, edges_sample_size, False)
        obs_sample = obs_edges.loc[altnodes_sample].reset_index(level=0, drop=True)
        exp_sample = freqs_nodes.loc[obs_sample.index]
        
        obs_sample_cnt = obs_sample.sum()
        exp_sample_cnt = exp_sample.sum()

        assert (obs_sample_cnt.index == exp_sample_cnt.index).all()

        sample_spectra = obs_sample_cnt / exp_sample_cnt
        spectra.append(sample_spectra)

    return pd.DataFrame(spectra).fillna(0.)


def calc_edgewise_spectra(
        obs: pd.DataFrame, exp: pd.DataFrame, 
        nmtypes_cutoff=10, nobs_cutoff=None, 
        collapse_to_12=False, scale=True, 
        both_12_and_192=False,
        nobs_cuttof=10,
    ):
    """
    Calculate per-branch (edge-wise) mutational spectra.

    For each tree branch the observed mutation counts are divided by the
    expected frequencies of the reference (parent) node, yielding a
    branch-specific spectrum.

    Arguments
    ---------
    obs: pd.DataFrame
        Observed mutations.  Either a pre-pivoted wide DataFrame with 192
        SBS columns and a ``(RefNode, AltNode)`` MultiIndex, or a long-format
        DataFrame with columns ``'RefNode'``, ``'AltNode'``, ``'Mut'``, and
        ``'ProbaFull'``.
    exp: pd.DataFrame
        Expected mutation frequencies.  Either a pre-pivoted wide DataFrame
        with 192 SBS columns and a ``Node`` index, or a long-format DataFrame
        with columns ``'Node'``, ``'Mut'``, and ``'Proba'``.
    nmtypes_cutoff: int
        Minimum number of distinct mutation types a branch must have to be
        retained (only applied when ``collapse_to_12=False``).
    nobs_cutoff: int
        Minimum total observed mutations a branch must have to be retained
        (only applied when ``collapse_to_12=False``).
    nobs_cuttof: int
        Deprecated alias of ``nobs_cutoff`` kept for backward compatibility.
    collapse_to_12: bool
        If ``True`` collapse the 192-component spectra to 12 components before
        returning.
    scale: bool
        If ``True`` normalise each branch spectrum to sum to 1.
    both_12_and_192: bool
        If ``True`` return both 12- and 192-component spectra as a tuple
        ``(spectra12, spectra192)``.

    Return
    ------
    spectra: pd.DataFrame or tuple[pd.DataFrame, pd.DataFrame]
        Branch-wise spectrum DataFrame (or tuple of two DataFrames when
        ``both_12_and_192=True``).
    """
    if nobs_cutoff is None:
        nobs_cutoff = nobs_cuttof
    if len(obs.columns) == 192 and \
            (obs.columns == possible_sbs192).all() and \
                (exp.columns == possible_sbs192).all():
        assert obs.index.names == ["RefNode", "AltNode"]
        assert exp.index.names == ["Node"]
        obs_edges = obs.copy()
        freqs_nodes = exp.copy()
        freqs_nodes.index.name = "RefNode"
    else:
        obs_edges = obs.groupby(["RefNode", "AltNode", "Mut"]).ProbaFull.sum().unstack()
        obs_edges = complete_sbs192_columns(obs_edges)
        freqs_nodes = exp.groupby(["Node", "Mut"]).Proba.sum().unstack()
        freqs_nodes.index.name = "RefNode"
        freqs_nodes = complete_sbs192_columns(freqs_nodes)

    if not collapse_to_12:
        obs_edges = obs_edges[((obs_edges > 0).sum(axis=1) >= nmtypes_cutoff) & \
                               (obs_edges.sum(axis=1) >= nobs_cutoff)]
    
    edges_df = obs_edges.index.to_frame(False)

    freqs_edges = edges_df.merge(freqs_nodes, on="RefNode")\
        .set_index(["RefNode", "AltNode"])[possible_sbs192]

    # some indexes can be deleted from freqs, so we must delete them from obs
    obs_edges = obs_edges.loc[freqs_edges.index]

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


def get_cossim(a: pd.DataFrame, b: pd.DataFrame):
    """
    Compute row-wise cosine similarity between two aligned DataFrames.

    Only rows present in both DataFrames (intersection of indices) are used.

    Arguments
    ---------
    a: pd.DataFrame
        First DataFrame; columns must match those of *b*.
    b: pd.DataFrame
        Second DataFrame; columns must match those of *a*.

    Return
    ------
    cossim: pd.Series
        Cosine similarity for each shared index, ranging from -1 to 1.
        Returns an empty Series if the indices do not overlap.
    """
    
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
    """
    Compute row-wise Euclidean distance between two aligned DataFrames.

    Only rows present in both DataFrames (intersection of indices) are used.

    Arguments
    ---------
    a: pd.DataFrame
        First DataFrame; columns must match those of *b*.
    b: pd.DataFrame
        Second DataFrame; columns must match those of *a*.

    Return
    ------
    d: pd.Series
        Euclidean distance for each shared index.
        Returns an empty Series if the indices do not overlap.
    """
    
    common_index = a.index.intersection(b.index)
    if len(common_index) == 0:
        return pd.Series()
    
    a = a.loc[common_index]
    b = b.loc[common_index]

    d = ((a - b) ** 2).sum(axis=1) ** 0.5
    return d
