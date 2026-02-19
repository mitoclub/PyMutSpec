"""
Functionality to plot mutational spectrums
"""

from typing import Iterable

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ordered_sbs12 = ["C>A", "G>T", "C>G", "G>C", "C>T", "G>A", 
                 "T>A", "A>T", "T>C", "A>G", "T>G", "A>C"]
ordered_sbs192 = [
    'A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 
    'C[C>A]G', 'C[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 
    'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[G>T]T', 'G[G>T]T', 
    'C[G>T]T', 'A[G>T]T', 'T[G>T]G', 'G[G>T]G', 'C[G>T]G', 'A[G>T]G', 
    'T[G>T]C', 'G[G>T]C', 'C[G>T]C', 'A[G>T]C', 'T[G>T]A', 'G[G>T]A', 
    'C[G>T]A', 'A[G>T]A', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 
    'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'G[C>G]A', 'G[C>G]C',
    'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 
    'T[G>C]T', 'G[G>C]T', 'C[G>C]T', 'A[G>C]T', 'T[G>C]G', 'G[G>C]G', 
    'C[G>C]G', 'A[G>C]G', 'T[G>C]C', 'G[G>C]C', 'C[G>C]C', 'A[G>C]C', 
    'T[G>C]A', 'G[G>C]A', 'C[G>C]A', 'A[G>C]A', 'A[C>T]A', 'A[C>T]C', 
    'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 
    'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 
    'T[C>T]G', 'T[C>T]T', 'T[G>A]T', 'G[G>A]T', 'C[G>A]T', 'A[G>A]T', 
    'T[G>A]G', 'G[G>A]G', 'C[G>A]G', 'A[G>A]G', 'T[G>A]C', 'G[G>A]C', 
    'C[G>A]C', 'A[G>A]C', 'T[G>A]A', 'G[G>A]A', 'C[G>A]A', 'A[G>A]A', 
    'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 
    'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 
    'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[A>T]T', 'G[A>T]T', 
    'C[A>T]T', 'A[A>T]T', 'T[A>T]G', 'G[A>T]G', 'C[A>T]G', 'A[A>T]G', 
    'T[A>T]C', 'G[A>T]C', 'C[A>T]C', 'A[A>T]C', 'T[A>T]A', 'G[A>T]A', 
    'C[A>T]A', 'A[A>T]A', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 
    'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 
    'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 
    'T[A>G]T', 'G[A>G]T', 'C[A>G]T', 'A[A>G]T', 'T[A>G]G', 'G[A>G]G', 
    'C[A>G]G', 'A[A>G]G', 'T[A>G]C', 'G[A>G]C', 'C[A>G]C', 'A[A>G]C', 
    'T[A>G]A', 'G[A>G]A', 'C[A>G]A', 'A[A>G]A', 'A[T>G]A', 'A[T>G]C', 
    'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 
    'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 
    'T[T>G]G', 'T[T>G]T', 'T[A>C]T', 'G[A>C]T', 'C[A>C]T', 'A[A>C]T', 
    'T[A>C]G', 'G[A>C]G', 'C[A>C]G', 'A[A>C]G', 'T[A>C]C', 'G[A>C]C', 
    'C[A>C]C', 'A[A>C]C', 'T[A>C]A', 'G[A>C]A', 'C[A>C]A', 'A[A>C]A'
]
ordered_sbs192_kp = ordered_sbs192

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
color_mapping192 = {}
for _sbs192 in ordered_sbs192:
    _sbs12 = _sbs192[2:5]
    color_mapping192[_sbs192] = color_mapping12[_sbs12]


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


def plot_mutspec(
        mutspec: pd.DataFrame,
        spectra_col="MutSpec",
        title="Spectrum",
        ylabel=None,
        figsize=None,
        style="bar",
        sbs_kind=12,
        sbs_order=None,
        savepath=None,
        fontname=None,
        ticksize=8,
        titlesize=14,
        ylabelsize=12,
        bar_width=0.8,
        show=True,
        dpi=300,
        ax=None,
        **kwargs,
    ):
    """
    General plotting function for mutational spectra supporting 
    12- and 192-component spectra.
    """
    if "filepath" in kwargs:
        savepath = kwargs.pop("filepath")

    is_192 = (sbs_kind == 192)

    # Defaults
    if figsize is None:
        figsize = (24, 8) if is_192 else (6, 4)

    if is_192:
        sbs_order = sbs_order or ordered_sbs192
        order = _prepare_nice_labels(sbs_order)
        palette = color_mapping192
        tick_rotation = 90
    else:
        order = sbs_order or ordered_sbs12
        palette = color_mapping12
        tick_rotation = 0

    ms = mutspec.copy()

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()

    if style == "bar":
        _cols = set(ms.columns)
        if 'MutSpec_median' in _cols and 'MutSpec_q05' in _cols and 'MutSpec_q95' in _cols:
            sns.barplot(
                ms, x='Mut', y='MutSpec', hue='Mut', legend=False, width=bar_width,
                order=order, palette=palette, err_kws={'linewidth': 1}, ax=ax, **kwargs)

            mutspec_index = ms.set_index('Mut')

            # build y and yerr aligned with order; support '_' separators in 192 order
            ymed, yerr_min, yerr_max = [], [], []
            for mt in order:
                if isinstance(mt, str) and mt.startswith('_'):
                    ymed.append(0.)
                    yerr_min.append(0.)
                    yerr_max.append(0.)
                else:
                    row = mutspec_index.loc[mt]
                    ymed.append(row['MutSpec_median'])
                    # convert to relative error for plotting
                    yerr_min.append(row['MutSpec_median'] - row['MutSpec_q05'])
                    yerr_max.append(row['MutSpec_q95'] - row['MutSpec_median'])

            yerr = [yerr_min, yerr_max]
            x = list(range(len(order)))
            ax.errorbar(x, ymed, yerr=yerr, fmt=".", color="gray", elinewidth=0.7,
                        capsize=(2 if is_192 else 4))
        else:
            sns.barplot(
                ms, x='Mut', y=spectra_col, hue='Mut', legend=False, width=bar_width,
                order=order, palette=palette, err_kws={'linewidth': 1}, ax=ax, **kwargs)

    elif style == "box":
        sns.boxplot(
            ms, x='Mut', y=spectra_col, hue='Mut', legend=False,
            order=order, palette=palette, ax=ax, **kwargs,
        )
    else:
        raise NotImplementedError

    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_ylabel(ylabel if ylabel else "", fontsize=ylabelsize, fontname=fontname)
    ax.set_xlabel("")
    ax.set_xlim(-0.5, len(order) - 0.5)

    if is_192:
        order_styled = ["" if x[0] == "_" else x for x in order]
        ax.set_xticks(range(len(order_styled)))
        ax.set_xticklabels(order_styled, rotation=tick_rotation, fontsize=ticksize, fontname=fontname)
    else:
        plt.xticks(fontsize=ticksize, fontname=fontname)

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
    return ax


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
        dpi=300,
        ax=None,
        **kwargs,
    ):
    # delegate to generic plotter for 12-component spectra
    return plot_mutspec(
        mutspec=mutspec,
        spectra_col=spectra_col,
        title=title,
        ylabel=ylabel,
        figsize=figsize,
        style=style,
        sbs_kind=12,
        sbs_order=ordered_sbs12,
        savepath=savepath,
        fontname=fontname,
        ticksize=ticksize,
        titlesize=titlesize,
        ylabelsize=ylabelsize,
        show=show,
        dpi=dpi,
        ax=ax,
        **kwargs,
    )


def plot_mutspec192(
        mutspec192: pd.DataFrame, 
        spectra_col="MutSpec", 
        title="Mutational spectrum", 
        ylabel=None, 
        figsize=(24, 8), 
        style="bar",
        sbs_order=ordered_sbs192_kp,
        savepath=None,
        fontname=None,
        ticksize=6,
        titlesize=16,
        ylabelsize=16,
        bar_width=0.6,
        show=True,
        dpi=300,
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
    # delegate to generic plotter for 192-component spectra
    return plot_mutspec(
        mutspec=mutspec192,
        spectra_col=spectra_col,
        title=title,
        ylabel=ylabel,
        figsize=figsize,
        style=style,
        sbs_kind=192,
        sbs_order=sbs_order,
        savepath=savepath,
        fontname=fontname,
        ticksize=ticksize,
        titlesize=titlesize,
        ylabelsize=ylabelsize,
        bar_width=bar_width,
        show=show,
        dpi=dpi,
        ax=ax,
        **kwargs,
    )


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
