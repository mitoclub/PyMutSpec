"""
Functionality to plot mutational spectrums
"""

from typing import Iterable

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..constants import possible_sbs192

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

# KK-style ordering: substitutions are grouped by base type and sorted using
# the SBS itself for kk_lbls, or its reverse-complement otherwise.
_kk_lbl_set = set("A>C A>G A>T C>T G>C G>T".split())
_transcriptor = str.maketrans("ACGT", "TGCA")


def _sbs192_rev_comp(sbs: str) -> str:
    """Return the reverse complement of a 192-component SBS string.

    The input must be a 7-character string of the form ``X[N>M]Y`` where
    ``X`` and ``Y`` are single flanking nucleotides and ``N>M`` is the
    substitution (e.g. ``'A[C>A]T'``).
    """
    return (sbs[-1] + sbs[1:-1] + sbs[0]).translate(_transcriptor)


ordered_sbs192_kk = sorted(
    possible_sbs192,
    key=lambda sbs: (sbs[2:5], sbs if sbs[2:5] in _kk_lbl_set else _sbs192_rev_comp(sbs)),
)

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
    """
    Build a display-friendly label list for 192-component SBS axes.

    A short separator string (underscores) is inserted between groups of
    substitutions that share the same base substitution type.

    Arguments
    ---------
    sbs192: iterable of str
        Ordered list of 192-component SBS strings (e.g. ``'A[C>A]C'``).
    kk: bool
        If ``True``, use the compact KK-style label format
        ``'CA: ACA'`` instead of the full COSMIC string.

    Return
    ------
    labels: list of str
        Label strings suitable for use as tick labels, with separator
        entries inserted between substitution groups.
    """
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
        labels_style="cosmic",
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

    Arguments
    ---------
    mutspec: pd.DataFrame
        Table containing at least a ``'Mut'`` column and a column named
        *spectra_col* with spectrum values.
    spectra_col: str
        Column name used for the y-axis values.
    title: str
        Plot title.
    ylabel: str or None
        Y-axis label.
    figsize: tuple or None
        Figure size ``(width, height)`` in inches.  Defaults to ``(24, 8)``
        for 192-component and ``(6, 4)`` for 12-component spectra.
    style: str
        ``'bar'`` for a bar plot or ``'box'`` for a box plot.
    sbs_kind: int
        Number of SBS components: ``12`` or ``192``.
    sbs_order: list or None
        Custom order for SBS labels on the x-axis.
    labels_style: str
        Style for 192-component tick labels.  One of:

        - ``'cosmic'`` – full COSMIC-format strings, e.g. ``'A[C>A]T'``
          (default)
        - ``'long'``   – same as ``'cosmic'``
        - ``'kk'``     – compact KK-style labels, e.g. ``'CA: ACT'``
    savepath: str or None
        File path to save the figure.  If ``None`` the figure is not saved.
    fontname: str or None
        Font family for all text elements.
    ticksize: int
        Font size for tick labels.
    titlesize: int
        Font size for the title.
    ylabelsize: int
        Font size for the y-axis label.
    bar_width: float
        Width of bars in bar plots.
    show: bool
        If ``True`` call ``plt.show()``; otherwise close the figure.
    dpi: int
        Resolution for saved figures.
    ax: matplotlib.axes.Axes or None
        Axes to draw on.  A new figure is created when ``None``.
    **kwargs
        Additional keyword arguments forwarded to the underlying seaborn
        plot function.

    Return
    ------
    ax: matplotlib.axes.Axes
        The axes containing the plot.
    """
    if "filepath" in kwargs:
        savepath = kwargs.pop("filepath")

    is_192 = (sbs_kind == 192)
    kk_labels = (labels_style == "kk")

    # Defaults
    if figsize is None:
        figsize = (24, 8) if is_192 else (6, 4)

    ms = mutspec.copy()

    if is_192:
        sbs_order = sbs_order or ordered_sbs192
        order = _prepare_nice_labels(sbs_order, kk=kk_labels)
        if kk_labels:
            # Map COSMIC SBS strings to KK-style display labels and rebuild palette
            sbs_to_kk = {
                sbs: sbs[2] + sbs[4] + ": " + sbs[0] + sbs[2] + sbs[-1]
                for sbs in possible_sbs192
            }
            ms['Mut'] = ms['Mut'].map(sbs_to_kk)
            palette = {sbs_to_kk[sbs]: color_mapping192[sbs] for sbs in possible_sbs192}
        else:
            palette = color_mapping192
        tick_rotation = 90
    else:
        order = sbs_order or ordered_sbs12
        palette = color_mapping12
        tick_rotation = 0

    created_fig = False
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()
        created_fig = True
    else:
        fig = ax.figure

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
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, fontsize=ticksize, fontname=fontname)

    if savepath is not None:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    elif created_fig:
        plt.close(fig)
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
    """
    Plot a 12-component mutational spectrum.

    A convenience wrapper around :func:`plot_mutspec` for 12-component
    (SBS12) spectra.

    Arguments
    ---------
    mutspec: pd.DataFrame
        Table containing at least a ``'Mut'`` column with 12-component SBS
        codes and a column named *spectra_col* with spectrum values.
    spectra_col: str
        Column name used for the y-axis values.
    title: str
        Plot title.
    ylabel: str or None
        Y-axis label.  Defaults to an empty string when ``None``.
    figsize: tuple
        Figure size ``(width, height)`` in inches.
    style: str
        ``'bar'`` for a bar plot or ``'box'`` for a box plot.
    savepath: str or None
        File path to save the figure.  If ``None`` the figure is not saved.
    fontname: str or None
        Font family for all text elements.
    ticksize: int
        Font size for tick labels.
    titlesize: int
        Font size for the title.
    ylabelsize: int
        Font size for the y-axis label.
    show: bool
        If ``True`` call ``plt.show()``; otherwise close the figure.
    dpi: int
        Resolution for saved figures.
    ax: matplotlib.axes.Axes or None
        Axes to draw on.  A new figure is created when ``None``.
    **kwargs
        Additional keyword arguments forwarded to the underlying seaborn
        plot function.

    Return
    ------
    ax: matplotlib.axes.Axes
        The axes containing the plot.
    """
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
        labels_style="cosmic",
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
    Plot a barplot of a 192-component mutational spectrum.

    Arguments
    ---------
    mutspec192: pd.DataFrame
        Table containing 192-component mutational spectrum for one or many
        species; all substitution types must be present in the table.
    spectra_col: str
        Column name used for the y-axis values.
    title: str
        Title on the plot.
    ylabel: str or None
        Y-axis label.
    figsize: tuple
        Figure size ``(width, height)`` in inches.
    style: str
        ``'bar'`` for a bar plot or ``'box'`` for a box plot.
    sbs_order: list or None
        Custom ordering of SBS192 labels on the x-axis.  Defaults to
        COSMIC ordering.
    labels_style: str
        Style for tick labels.  One of ``'cosmic'``/``'long'`` (full COSMIC
        strings) or ``'kk'`` (compact KK-style labels).
    savepath: str or None
        Path to output plot file.  If ``None`` no image is written.
    fontname: str or None
        Font family for all text elements.
    ticksize: int
        Font size for tick labels.
    titlesize: int
        Font size for the title.
    ylabelsize: int
        Font size for the y-axis label.
    bar_width: float
        Width of bars.
    show: bool
        If ``True`` call ``plt.show()``; otherwise close the figure.
    dpi: int
        Resolution for saved figures.
    ax: matplotlib.axes.Axes or None
        Axes to draw on.  A new figure is created when ``None``.
    **kwargs
        Additional keyword arguments forwarded to the underlying seaborn
        plot function.

    Return
    ------
    ax: matplotlib.axes.Axes
        The axes containing the plot.
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
        labels_style=labels_style,
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
