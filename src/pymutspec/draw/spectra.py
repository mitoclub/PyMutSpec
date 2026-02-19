"""
Functionality to plot mutational spectrums
"""

from typing import Iterable

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .sbs_orders import ordered_sbs192_kp

color_mapping12 = {
    "C>A": "deepskyblue",
    "G>T": "deepskyblue",
    "C>G": "black",
    "G>C": "black",
    "C>T": "red",
    "G>A": "red",
    "T>A": "silver",
    "A>T": "silver",
    "T>C": "yellowgreen",
    "A>G": "yellowgreen",
    "T>G": "pink",
    "A>C": "pink",
}
sbs12_ordered = ["C>A", "G>T", "C>G", "G>C", "C>T", "G>A", 
                 "T>A", "A>T", "T>C", "A>G", "T>G", "A>C"]

color_mapping192 = {}
for sbs, clr in color_mapping12.items():
    for nuc1 in "ACGT":
        for nuc2 in "ACGT":
            sbs192 = f'{nuc1}[{sbs}]{nuc2}'
            color_mapping192[sbs192] = clr


def plot_mutspec12(
        mutspec: pd.DataFrame, 
        spectra_col="MutSpec", 
        title="Spectrum", 
        ylabel=None, 
        figsize=(6, 4), 
        style="bar", 
        savepath=None, 
        fontname=None,
        ticksize=8,
        titlesize=14,
        ylabelsize=12,
        show=True,
        ax=None,
        **kwargs,
    ):
    # TODO add checks of mutspec12
    # TODO add description to all plot* functions

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()

    if style == "bar":
        _cols = set(mutspec.columns)
        if 'MutSpec_median' in _cols and 'MutSpec_q05' in _cols and 'MutSpec_q95' in _cols:
            ax = sns.barplot(
                data=mutspec, x="Mut", y='MutSpec', hue='Mut', legend=False,
                order=sbs12_ordered, palette=color_mapping12, err_kws={'linewidth': 1},
                ax=ax, **kwargs,
            )
            
            mutspec_ordered = mutspec.set_index('Mut').loc[sbs12_ordered]
            mutspec_ordered['MutSpec_q05'] = mutspec_ordered['MutSpec_median'] - mutspec_ordered['MutSpec_q05']
            mutspec_ordered['MutSpec_q95'] -= mutspec_ordered['MutSpec_median']
            ax.errorbar(mutspec_ordered.index, mutspec_ordered['MutSpec_median'], 
                        yerr=mutspec_ordered[['MutSpec_q05', 'MutSpec_q95']].values.T, 
                        fmt=".", color="gray", elinewidth=0.7, capsize=4)
        else:
            ax = sns.barplot(
                data=mutspec, x="Mut", y=spectra_col, hue='Mut', legend=False,
                order=sbs12_ordered, palette=color_mapping12,
                ax=ax, **kwargs,
            )

    elif style == "box":
        ax = sns.boxplot(
                data=mutspec, x="Mut", y=spectra_col, hue='Mut', legend=False,
                order=sbs12_ordered, palette=color_mapping12,
                ax=ax, **kwargs,
            )
    else:
        raise NotImplementedError
    
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_ylabel(ylabel if ylabel else "", fontsize=ylabelsize, fontname=fontname)
    ax.set_xlabel("")
    plt.xticks(fontsize=ticksize, fontname=fontname)

    if savepath is not None:
        plt.savefig(savepath, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
    return ax


def _prepare_nice_labels(sbs192: Iterable[str], kk=False):
    _nice_sbs = []
    prev = None
    for i, sbs in enumerate(sbs192, 1):
        if prev is not None and sbs[2:5] != prev[2:5]:
            _nice_sbs.append("_" * (i // 10))
        sbs_nice = sbs[2] + sbs[4] + ": " + sbs[0] + sbs[2] + sbs[-1] if kk else sbs
        _nice_sbs.append(sbs_nice)
        prev = sbs
    return _nice_sbs


def plot_mutspec192(
        mutspec192: pd.DataFrame, 
        spectra_col="MutSpec", 
        title="Mutational spectrum", 
        ylabel=None, 
        figsize=(24, 8), 
        style="bar",
        sbs_order='full_192',
        savepath=None,
        fontname=None,
        ticksize=6,
        titlesize=16,
        ylabelsize=16,
        bar_width=0.6,
        show=True,
        ax=None,
        **kwargs,
    ):
    """
    Plot barblot of given mutational spectrum calculated from single nucleotide substitutions

    Arguments
    ---------
    mutspec192: pd.DataFrame
        table, containing 192 component mutational spectrum for one or many species, all substitutions must be presented in the table
    title: str, default = 'Mutational spectrum'
        Title on the plot
    savepath: str, default = None
        Path to output plot file. If None no images will be written
    """
    if "filepath" in kwargs:
        savepath = kwargs["filepath"]
        kwargs.pop("filepath")

    if sbs_order == 'full_192':
        sbs_order = ordered_sbs192_kp

    # TODO add checks of mutspec192
    ms192 = mutspec192.copy()
    order = _prepare_nice_labels(sbs_order, kk=False)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()
    
    if style == "bar":
        _cols = set(ms192.columns)
        if 'MutSpec_median' in _cols and 'MutSpec_q05' in _cols and 'MutSpec_q95' in _cols:
            sns.barplot(
                ms192, x='Mut', y='MutSpec', hue='Mut', legend=False, width=bar_width,
                order=order, palette=color_mapping192, err_kws={'linewidth': 1}, ax=ax, **kwargs)
            
            mutspec_ordered = ms192.set_index('Mut')
            mutspec_ordered['MutSpec_q05'] = mutspec_ordered['MutSpec_median'] - mutspec_ordered['MutSpec_q05']
            mutspec_ordered['MutSpec_q95'] -= mutspec_ordered['MutSpec_median']

            ymed, yerr_min, yerr_max = [], [], []
            for mt in order:
                if mt[0] == '_':
                    ymed.append(0.)
                    yerr_min.append(0.)
                    yerr_max.append(0.)
                else:
                    ymed.append(mutspec_ordered.loc[mt, 'MutSpec_median'])
                    yerr_min.append(mutspec_ordered.loc[mt, 'MutSpec_q05'])
                    yerr_max.append(mutspec_ordered.loc[mt, 'MutSpec_q95'])
            yerr = [yerr_min, yerr_max]
            x = list(range(len(order)))
            ax.errorbar(x, ymed, yerr=yerr, fmt=".", color="gray", elinewidth=0.7, capsize=2)
        else:
            sns.barplot(
                ms192, x='Mut', y=spectra_col, hue='Mut', legend=False, width=bar_width,
                order=order, palette=color_mapping192, err_kws={'linewidth': 1}, ax=ax, **kwargs)
    elif style == "box":
        sns.boxplot(
            ms192, x='Mut', y=spectra_col, hue='Mut', legend=False, 
            order=order, palette=color_mapping192, ax=ax, **kwargs,
        )
    
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_ylabel(ylabel if ylabel else "", fontsize=ylabelsize, fontname=fontname)
    ax.set_xlabel("")

    order_styled = ["" if x[0] == "_" else x for x in order]
    ax.set_xticks(range(len(order_styled)))
    ax.set_xticklabels(order_styled, rotation=90, fontsize=ticksize, fontname=fontname)

    plt.xlim(-0.5, len(order) - 0.5)

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
    return ax


# def plot_mutspec192kk(mutspec192: pd.DataFrame, ylabel="MutSpec", title="Mutational spectrum", show=True, figsize=(24, 6), filepath=None):
    # from .sbs_orders import ordered_sbs192_kk
    # ms192 = mutspec192.copy()
#     ms192["long_lbl"] = ms192.Mut.str.get(2) + ms192.Mut.str.get(4) + ": " + ms192.Mut.str.get(0) + ms192.Mut.str.get(2) + ms192.Mut.str.get(-1)
#     fig = plt.figure(figsize=figsize)
#     ax = fig.add_subplot(111)
#     ax.grid(axis="y", alpha=.7, linewidth=0.5)
#     order = _prepare_nice_labels(ordered_sbs192_kk, True)
#     sns.barplot(
#         x="long_lbl", y=ylabel, data=ms192,
#         order=order, 
#         errwidth=1, ax=fig.gca(), 
#     )
#     plt.xticks(rotation=90, fontsize=7, fontname=None)
#     ax.set_title(title)
#     ax.set_xlabel("")
#     ax.set_ylabel("Mutational spectrum")

#     def _coloring192kk():
#         colors = "red yellow lime blue".split()
#         while True:
#             for clr in colors:
#                 yield clr

#     # map colors to bars
#     clrs_iterator = _coloring192kk()
#     for bar, sbs in zip(ax.patches, order):
#         if len(sbs):
#             bar.set_color(next(clrs_iterator))
#             bar.set_alpha(alpha=0.9)
#         bar.set_width(0.3)
#     if filepath is not None:
#         plt.savefig(filepath, dpi=300, bbox_inches="tight")
#     if show:
#         plt.show()
#     else:
#         plt.close()
