"""
Functionality to plot mutational spectrums
"""

from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pymutspec.constants import possible_sbs192, possible_sbs96
from .sbs_orders import ordered_sbs192_kk, ordered_sbs192_kp

color_mapping6 = {
    "C>A": "deepskyblue",
    "C>G": "black",
    "C>T": "red",
    "T>A": "silver",
    "T>C": "green",
    "T>G": "pink",
}
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
sbs12_ordered = ["C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C"]
_smpl = pd.DataFrame({"Mut": possible_sbs192})
_smpl["MutBase"] = _smpl.Mut.str.slice(2, 5)
_smpl["Context"] = _smpl.Mut.str.get(0) + _smpl.Mut.str.get(2) + _smpl.Mut.str.get(-1)
colors192 = _smpl.sort_values(["MutBase", "Context"])["MutBase"].map(color_mapping12).values
colors12 = [color_mapping12[sbs] for sbs in sbs12_ordered]


def _prepare_nice_labels(sbs192: Iterable[str], kk=False):
    _nice_sbs = []
    prev = None
    for sbs in sbs192:
        if prev is not None and sbs[2:5] != prev[2:5]:
            _nice_sbs.append("")
        sbs_nice = sbs[2] + sbs[4] + ": " + sbs[0] + sbs[2] + sbs[-1] if kk else sbs
        _nice_sbs.append(sbs_nice)
        prev = sbs
    return _nice_sbs


def _coloring192kk():
    colors = "red yellow lime blue".split()
    while True:
        for clr in colors:
            yield clr


def plot_mutspec12(mutspec: pd.DataFrame, ylabel="MutSpec", title="Full mutational spectrum", figsize=(6, 4), style="bar", show=True, savepath=None, **kwargs):
    # TODO add checks of mutspec12
    # TODO add description to all plot* functions
    if style == "bar":
        plot_func = sns.barplot
    elif style == "box":
        plot_func = sns.boxplot
    else:
        raise NotImplementedError

    fig = plt.figure(figsize=figsize)
    ax = plot_func(x="Mut", y=ylabel, data=mutspec, order=sbs12_ordered, ax=fig.gca(), **kwargs)
    ax.grid(axis="y", alpha=.7, linewidth=0.5)

    # map colors to bars
    for bar, clr in zip(ax.patches, colors12):
        bar.set_color(clr)

    ax.set_title(title)
    ax.set_ylabel("")
    ax.set_xlabel("")

    if savepath is not None:
        plt.savefig(savepath)
    if show:
        plt.show()
    else:
        plt.close()
    return ax


def plot_mutspec192(
        mutspec192: pd.DataFrame, 
        ylabel="MutSpec", 
        title="Mutational spectrum", 
        figsize=(24, 8), 
        style="bar",
        labels_style="cosmic",
        sbs_order=ordered_sbs192_kp,
        savepath=None,
        fontsize=6,
        titlesize=16,
        fontname="Times New Roman",
        show=True, 
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
    labels_style: str, default = 'cosmic'
        'cosmic': A[C>T]T, 'long': CT: ACT
    """
    if "filepath" in kwargs:
        savepath = kwargs["filepath"]
        print("savepath =", savepath)
        kwargs.pop("filepath")
    
    # TODO add checks of mutspec192
    ms192 = mutspec192.copy()
    if labels_style == "long":
        ms192["long_lbl"] = ms192.Mut.str.get(2) + ms192.Mut.str.get(4) + ": " + ms192.Mut.str.get(0) + ms192.Mut.str.get(2) + ms192.Mut.str.get(-1)
        order = _prepare_nice_labels(sbs_order, kk=True)
        x_col = "long_lbl"
    elif labels_style == "cosmic":
        order = _prepare_nice_labels(sbs_order, kk=False)
        x_col = "Mut"
    else:
        raise ValueError("Available labels_style are: 'cosmic' and 'long'")
    fig = plt.figure(figsize=figsize)
    if style == "bar":
        ax = sns.barplot(
            x=x_col, y=ylabel, data=ms192, order=order, errwidth=1, ax=fig.gca(), **kwargs,
        )
    elif style == "box":
        ax = sns.boxplot(
            x=x_col, y=ylabel, data=ms192, order=order, ax=fig.gca(), **kwargs,
        )
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_xlabel("")
    ax.set_ylabel("")
    # map colors to bars
    width = 0.4
    shift = None
    for bar, sbs in zip(ax.patches, order):
        if len(sbs):
            s = sbs[0] + ">" + sbs[1] if labels_style == "long" else sbs[2:5]
            bar.set_color(color_mapping12[s])
            bar.set_alpha(alpha=0.9)
        if style == "bar":
            if not shift:
                # calculate one time instead of 192
                shift = (bar.get_width() - width) / 2
            bar.set_width(width)
            bar.set_x(bar.get_x() + shift)

    plt.xticks(rotation=90, fontsize=fontsize, fontname=fontname)

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def plot_mutspec96(
        spectra: pd.DataFrame, 
        ylabel="MutSpec", 
        title="Mutational spectrum", 
        figsize=(15, 5), 
        style="bar",
        labels_style="cosmic",
        sbs_order=possible_sbs96,
        savepath=None,
        fontsize=9,
        titlesize=16,
        fontname="Times New Roman",
        show=True, 
        **kwargs,
    ):
    """
    Plot barblot of given mutational spectrum calculated from single nucleotide substitutions

    Arguments
    ---------
    spectra: pd.DataFrame
        table, containing 192 component mutational spectrum for one or many species, all substitutions must be presented in the table
    title: str, default = 'Mutational spectrum'
        Title on the plot
    savepath: str, default = None
        Path to output plot file. If None no images will be written
    labels_style: str, default = 'cosmic'
        'cosmic': A[C>T]T, 'long': CT: ACT
    """
    if "filepath" in kwargs:
        savepath = kwargs["filepath"]
        print("savepath =", savepath)
        kwargs.pop("filepath")
    
    # TODO add checks of spectra
    ms96 = spectra.copy()
    if labels_style == "long":
        ms96["long_lbl"] = ms96.Mut.str.get(2) + ms96.Mut.str.get(4) + ": " + ms96.Mut.str.get(0) + ms96.Mut.str.get(2) + ms96.Mut.str.get(-1)
        order = _prepare_nice_labels(sbs_order, kk=True)
        x_col = "long_lbl"
    elif labels_style == "cosmic":
        order = _prepare_nice_labels(sbs_order, kk=False)
        x_col = "Mut"
    else:
        raise ValueError("Available labels_style are: 'cosmic' and 'long'")
    fig = plt.figure(figsize=figsize)
    if style == "bar":
        ax = sns.barplot(
            x=x_col, y=ylabel, data=ms96, order=order, errwidth=1, ax=fig.gca(), **kwargs,
        )
    elif style == "box":
        ax = sns.boxplot(
            x=x_col, y=ylabel, data=ms96, order=order, ax=fig.gca(), **kwargs,
        )
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_xlabel("")
    ax.set_ylabel("")
    # map colors to bars
    width = 0.4
    shift = None
    for bar, sbs in zip(ax.patches, order):
        if len(sbs):
            s = sbs[0] + ">" + sbs[1] if labels_style == "long" else sbs[2:5]
            bar.set_color(color_mapping6[s])
            bar.set_alpha(alpha=0.9)
        if style == "bar":
            if not shift:
                # calculate one time instead of 192
                shift = (bar.get_width() - width) / 2
            bar.set_width(width)
            bar.set_x(bar.get_x() + shift)

    plt.xticks(rotation=90, fontsize=fontsize, fontname=fontname)

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def plot_mutspec192kk(mutspec192: pd.DataFrame, ylabel="MutSpec", title="Mutational spectrum", show=True, figsize=(24, 6), filepath=None):
    ms192 = mutspec192.copy()
    ms192["long_lbl"] = ms192.Mut.str.get(2) + ms192.Mut.str.get(4) + ": " + ms192.Mut.str.get(0) + ms192.Mut.str.get(2) + ms192.Mut.str.get(-1)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    order = _prepare_nice_labels(ordered_sbs192_kk, True)
    sns.barplot(
        x="long_lbl", y=ylabel, data=ms192,
        order=order, 
        errwidth=1, ax=fig.gca(), 
    )
    plt.xticks(rotation=90, fontsize=7, fontname="Times New Roman")
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Mutational spectrum")
    # map colors to bars
    clrs_iterator = _coloring192kk()
    for bar, sbs in zip(ax.patches, order):
        if len(sbs):
            bar.set_color(next(clrs_iterator))
            bar.set_alpha(alpha=0.9)
        bar.set_width(0.3)
    if filepath is not None:
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
