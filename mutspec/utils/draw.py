from itertools import groupby

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .constants import possible_sbs12, possible_sbs192


def plot_mutspec12(mutspec: pd.DataFrame, ylabel="MutSpec", title="Full mutational spectra", savepath=None):
    fig = plt.figure(figsize=(24, 12))
    ax = fig.add_subplot(111)
    ax = sns.barplot(x="Mut", y=ylabel, data=mutspec, order=possible_sbs12, ax=fig.gca())
    ax.set_title(title)
    if savepath is not None:
        plt.savefig(savepath)
    return ax


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
        'family': 'cursive',
        'color':  'black',
        'weight': 'normal',
        'size': 7,
    }
    rotation = 90
    ypos = -.05
    scale = 1./df.index.size
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


def plot_mutspec192(mutspec192: pd.DataFrame, title="Mutational spectra", filepath=None):
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
    mutspec192 = mutspec192.copy()
    mutspec192["MutBase"] = mutspec192.Mut.str.slice(2, 5)
    mutspec192["Context"] = mutspec192.Mut.str.get(0) + mutspec192.Mut.str.get(2) + mutspec192.Mut.str.get(-1)

    df = mutspec192.groupby(["MutBase", "Context"]).mean()
    fig = plt.figure(figsize=(24, 12))
    ax = fig.add_subplot(111)
    sns.barplot(x="Mut", y="MutSpec", data=mutspec192,
                order=possible_sbs192, errwidth=1, ax=fig.gca())

    labels = ['' for _ in ax.get_xticklabels()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('')
    ax.set_title(title)
    __label_group_bar_table(ax, df)
    fig.subplots_adjust(bottom=0.1 * df.index.nlevels)
    if filepath is not None:
        plt.savefig(filepath)
    plt.show()
