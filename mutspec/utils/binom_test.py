from typing import Dict, Union

import pandas as pd
import scipy.stats
import statsmodels.api as sm

reciprocal_pairs = [
    ("A>C", "C>A"),
    ("A>G", "G>A"),
    ("A>T", "T>A"),
    ("C>G", "G>C"),
    ("C>T", "T>C"),
    ("G>T", "T>G"),
]
complementary_pairs = [
    ("A>C", "T>G"),
    ("A>G", "T>C"),
    ("A>T", "T>A"),
    ("C>G", "G>C"),
    ("C>T", "G>A"),
    ("G>T", "C>A"),
]


def asterisk_for_pvalue(pval: float) -> str:
    if pval < 0.001:
        asterisk = "***"
    elif pval < 0.01:
        asterisk = "**"
    elif pval < 0.05:
        asterisk = "*"
    else:
        asterisk = ""
    return asterisk


def binom_testing(
        mut_counts: Dict[str, Union[int, float]], mode="complementary", 
        scale=False, adjust_pvals=True, adj_method="fdr_bh"
    ):
    """
    Apply binomial test to pairs of mutations

    Arguments
    ---------
    mut_counts: Dict[str, Union[int, float]]
        sbs counts (12 numbers)
    mode: string
        How to apply pairwise comparison to mutations: complementary (C>T and G>A ...) or reciprocal (C>T and T>C ...)
    scale: bool
        Scale counts (divide to lowest value in pair) in case of counts < 0.7 or not, default False
    adjust_pvals: bool
        Adjust pvalues or not, default True
    adj_method: str
        Method used for adjusting. See sm.stats.multipletests method documentation, default 'fdr_bh'

    Return
    ------
        pd.DataFrame with statistics values
    """
    if mode == "complementary":
        pairs = complementary_pairs
    elif mode == "reciprocal":
        pairs = reciprocal_pairs
    else:
        raise NotImplementedError()

    data = []
    for mut1, mut2 in pairs:
        n1, n2 = mut_counts[mut1], mut_counts[mut2]
        if n1 < n2:
            n1, n2 = n2, n1
            mut1, mut2 = mut2, mut1
        ratio = n1 / n2
        if scale and (n1 < 0.7 or n2 < 0.7):
            sc_factor = min(n1, n2)
            n1, n2 = n1 / sc_factor, n2 / sc_factor
        n1, n2 = round(n1), round(n2)
        res = scipy.stats.binomtest(n1, n1 + n2, p=0.5)
        pval = res.pvalue
        data.append((mut1, mut2, ratio, pval))
    cols = ["mut1", "mut2", "ratio", "pval"]
    df = pd.DataFrame(data, columns=cols)
    if adjust_pvals:
        _, qval, _, _ = sm.stats.multipletests(df["pval"].values, method="fdr_bh")  # adjust pval
        df["pval_adj"] = qval
        df["asterics"] = list(map(asterisk_for_pvalue, qval))
    else:
        df["asterics"] = list(map(asterisk_for_pvalue, df["pval"].values))
    return df
