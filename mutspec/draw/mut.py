"""
Functionality to plot mutational spectrums
"""

from itertools import groupby

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..constants import possible_sbs12, possible_sbs192, ordered_sbs192


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
    "G>A": "salmon",
    "T>A": "silver",
    "A>T": "silver",
    "T>C": "green",
    "A>G": "green",
    "T>G": "pink",
    "A>C": "pink",
}
_smpl = pd.DataFrame({"Mut": possible_sbs192})
_smpl["MutBase"] = _smpl.Mut.str.slice(2, 5)
_smpl["Context"] = _smpl.Mut.str.get(0) + _smpl.Mut.str.get(2) + _smpl.Mut.str.get(-1)
colors192 = _smpl.sort_values(["MutBase", "Context"])["MutBase"].map(color_mapping12).values
colors12 = [color_mapping12[sbs] for sbs in possible_sbs12]


def __prepare_nice_kk_order():
    _nice_order_kk = []
    prev = None
    for sbs in ordered_sbs192:
        if prev is not None and sbs[2:5] != prev[2:5]:
            _nice_order_kk.append("")
        _nice_order_kk.append(sbs[2] + sbs[4] + ": " + sbs[0] + sbs[2] + sbs[-1])
        prev = sbs
    return _nice_order_kk


_nice_order_kk = __prepare_nice_kk_order()


def _coloring192():
    colors = "red yellow lime blue".split()
    while True:
        for clr in colors:
            yield clr


def plot_mutspec12(mutspec: pd.DataFrame, ylabel="MutSpec", title="Full mutational spectra", savepath=None):
    # TODO add checks of mutspec12
    # TODO add description to all plot* functions
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax = sns.barplot(x="Mut", y=ylabel, data=mutspec, order=possible_sbs12, ax=fig.gca())

    # map colors to bars
    for bar, clr in zip(ax.patches, colors12):
        bar.set_color(clr)

    ax.set_title(title)
    if savepath is not None:
        plt.savefig(savepath)
    plt.show()


def __add_line(ax, xpos, ypos):
    line = plt.Line2D([xpos, xpos], [ypos + .1, ypos],
                      transform=ax.transAxes, color='black', linewidth=1)
    line.set_clip_on(False)
    ax.add_line(line)


def __label_len(my_index, level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k, g in groupby(labels)]


def __label_group_bar_table(ax, df):
    font = {
        'color':  'black',
        'weight': 'normal',
        'size': 7,
    }
    rotation = 90
    ypos = -.05
    scale = 1. / df.index.size
    for level in range(df.index.nlevels)[::-1]:
        if level == 0:
            rotation = 0
            font['size'] = 12

        pos = 0
        for label, rpos in __label_len(df.index, level):
            lxpos = (pos + .5 * rpos)*scale
            ax.text(lxpos, ypos, label, ha='center', rotation=rotation,
                    fontdict=font, transform=ax.transAxes)
            if level == 0:
                __add_line(ax, pos*scale, ypos)
            pos += rpos
        if level == 0:
            __add_line(ax, pos*scale, ypos)
        ypos -= .05


def plot_mutspec192(mutspec192: pd.DataFrame, ylabel="MutSpec", title="Mutational spectra", figsize=(20, 12), filepath=None):
    """
    Plot barblot of given mutational spectra calculated from single nucleotide substitutions

    Arguments
    ---------
    mutspec192: pd.DataFrame
        table, containing 192 component mutational spectra for one or many species, all substitutions must be presented in the table
    title: str, default = 'Mutational spectra'
        Title on the plot
    filepath: str, default = None
        Path to output plot file. If None no images will be written
    """
    # TODO add checks of mutspec192
    mutspec192 = mutspec192.copy()
    mutspec192["MutBase"] = mutspec192.Mut.str.slice(2, 5)
    mutspec192["Context"] = mutspec192.Mut.str.get(0) + mutspec192.Mut.str.get(2) + mutspec192.Mut.str.get(-1)

    df = mutspec192.groupby(["MutBase", "Context"]).mean()
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    sns.barplot(
        x="Mut", y=ylabel, data=mutspec192,
        order=ordered_sbs192, 
        errwidth=1, ax=fig.gca(),
    )
    # map colors to bars
    for bar, clr in zip(ax.patches, colors192):
        bar.set_color(clr)
        bar.set_width(0.3)

    labels = ['' for _ in ax.get_xticklabels()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('')
    ax.set_title(title)
    __label_group_bar_table(ax, df)
    fig.subplots_adjust(bottom=0.1 * df.index.nlevels)
    if filepath is not None:
        plt.savefig(filepath)
    plt.show()


def plot_mutspec192kk(mutspec192: pd.DataFrame, ylabel="MutSpec", title="Mutational spectra", figsize=(24, 6), filepath=None):
    ms192 = mutspec192.copy()
    ms192["long_lbl"] = ms192.Mut.str.get(2) + ms192.Mut.str.get(4) + ": " + ms192.Mut.str.get(0) + ms192.Mut.str.get(2) + ms192.Mut.str.get(-1)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    sns.barplot(
        x="long_lbl", y=ylabel, data=ms192,
        order=_nice_order_kk, errwidth=1, ax=fig.gca()
    )
    plt.xticks(rotation=90, fontsize=7)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Mutational spectrum")
    # map colors to bars
    clrs_iterator = _coloring192()
    for bar, sbs in zip(ax.patches, _nice_order_kk):
        if len(sbs):
            bar.set_color(next(clrs_iterator))
            bar.set_alpha(alpha=0.9)
        bar.set_width(0.3)
    if filepath is not None:
        plt.savefig(filepath)
    plt.show()
