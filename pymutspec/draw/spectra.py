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


def plot_mutspec12(
        mutspec: pd.DataFrame, 
        spectra_col="MutSpec", 
        title="Full mutational spectrum", 
        ylabel=None, 
        figsize=(6, 4), 
        style="bar", 
        savepath=None, 
        fontname=None,
        ticksize=8,
        titlesize=14,
        ylabelsize=12,
        show=True, 
        **kwargs,
    ):
    # TODO add checks of mutspec12
    # TODO add description to all plot* functions
    if style == "bar":
        plot_func = sns.barplot
    elif style == "box":
        plot_func = sns.boxplot
    else:
        raise NotImplementedError

    fig = plt.figure(figsize=figsize)
    ax = plot_func(x="Mut", y=spectra_col, data=mutspec, order=sbs12_ordered, ax=fig.gca(), **kwargs)
    ax.grid(axis="y", alpha=.7, linewidth=0.5)

    # map colors to bars
    for bar, clr in zip(ax.patches, colors12):
        bar.set_color(clr)

    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_ylabel(ylabel if ylabel else "", fontsize=ylabelsize, fontname=fontname)
    ax.set_xlabel("")

    plt.xticks(fontsize=ticksize, fontname=fontname)

    if savepath is not None:
        plt.savefig(savepath)
    if show:
        plt.show()
    else:
        plt.close()
    return ax


def plot_mutspec192(
        mutspec192: pd.DataFrame, 
        spectra_col="MutSpec", 
        title="Mutational spectrum", 
        ylabel=None, 
        figsize=(24, 8), 
        style="bar",
        labels_style="cosmic",
        sbs_order=ordered_sbs192_kp,
        savepath=None,
        fontname=None,
        ticksize=6,
        titlesize=16,
        ylabelsize=16,
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
            x=x_col, y=spectra_col, data=ms192, order=order, errwidth=1, ax=fig.gca(), **kwargs,
        )
    elif style == "box":
        ax = sns.boxplot(
            x=x_col, y=spectra_col, data=ms192, order=order, ax=fig.gca(), **kwargs,
        )
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    ax.set_title(title, fontsize=titlesize, fontname=fontname)
    ax.set_ylabel(ylabel if ylabel else "", fontsize=ylabelsize, fontname=fontname)
    ax.set_xlabel("")
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

    plt.xticks(rotation=90, fontsize=ticksize, fontname=fontname)

    if savepath is not None:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
    return ax


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
        fontname=None,
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
    return ax


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
    plt.xticks(rotation=90, fontsize=7, fontname=None)
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

def draw192mutspec(path_to_mut_data='../data/u_mutation_dists.filtered.csv',
                   save_path='../figures/MutSpecSyn192.pdf',
                   title='192 component mutation spectrum by \n synonymous mutations SARS-CoV-2',
                   parent_nucl_context='parent_nucl_context',
                   child_nucl_context='child_nucl_context',
                   aasub='syn',
                   replace_t_u=True,
                   add_n=True):
    """
    Function for drawing 192 component mutation spectrum
    Parameters:
        path_to_mut_data - path to your dataframe with mutations
        save_path - path to figure folder and figure name
        title - figure title
        parent_nucl_context - name of column from mutation dataset with mutation context(parent). Example - "taGgc"
        child_nucl_context - name of column from mutation dataset with mutation context(child). Example - "taTgc"
        aasub - If you want to specify your spectrum, you need to have column "AaSub" with mutation explanation - synonymous (S)/non synonymous (NS)/four-fold (FF)
            possible values = 'S'-syn+ff
                              'FF' - ff
                              'NS' - ns
        replace_t_u - If you want to replace T to U in your figure(True/False). Default = True
        add_n - If you need to add number of mutations in your figure(True/False). Default = True
    """

    df = pd.read_csv(path_to_mut_data)

    if aasub=='syn':
        df = df[df['AaSub'].isin(['S', 'FF'])]
    elif aasub=='ff':
        df = df[df['AaSub'].isin(['FF'])]
    elif aasub=='ns':
        df = df[df['AaSub'].isin(['NS'])]

    df['from'] = df[parent_nucl_context].astype(str).str[1:4].str.upper()
    df['to'] = df[child_nucl_context].astype(str).str[1:4].str.upper()
    df['from_to_context'] = df['from'] + '>' + df['to']

    dummy_list = ['AAA>ACA', 'AAA>AGA', 'AAA>ATA', 'AAC>ACC', 'AAC>AGC', 'AAC>ATC', 'AAG>ACG', 'AAG>AGG', 'AAG>ATG', 'AAT>ACT',
                   'AAT>AGT', 'AAT>ATT', 'ACA>AAA', 'ACA>AGA', 'ACA>ATA', 'ACC>AAC', 'ACC>AGC', 'ACC>ATC', 'ACG>AAG', 'ACG>AGG',
                   'ACG>ATG', 'ACT>AAT', 'ACT>AGT', 'ACT>ATT', 'AGA>AAA', 'AGA>ACA', 'AGA>ATA', 'AGC>AAC', 'AGC>ACC', 'AGC>ATC',
                   'AGG>AAG', 'AGG>ACG', 'AGG>ATG', 'AGT>AAT', 'AGT>ACT', 'AGT>ATT', 'ATA>AAA', 'ATA>ACA', 'ATA>AGA', 'ATC>AAC',
                   'ATC>ACC', 'ATC>AGC', 'ATG>AAG', 'ATG>ACG', 'ATG>AGG', 'ATT>AAT', 'ATT>ACT', 'ATT>AGT', 'CAA>CCA', 'CAA>CGA',
                   'CAA>CTA', 'CAC>CCC', 'CAC>CGC', 'CAC>CTC', 'CAG>CCG', 'CAG>CGG', 'CAG>CTG', 'CAT>CCT', 'CAT>CGT', 'CAT>CTT',
                   'CCA>CAA', 'CCA>CGA', 'CCA>CTA', 'CCC>CAC', 'CCC>CGC', 'CCC>CTC', 'CCG>CAG', 'CCG>CGG', 'CCG>CTG', 'CCT>CAT',
                   'CCT>CGT', 'CCT>CTT', 'CGA>CAA', 'CGA>CCA', 'CGA>CTA', 'CGC>CAC', 'CGC>CCC', 'CGC>CTC', 'CGG>CAG', 'CGG>CCG',
                   'CGG>CTG', 'CGT>CAT', 'CGT>CCT', 'CGT>CTT', 'CTA>CAA', 'CTA>CCA', 'CTA>CGA', 'CTC>CAC', 'CTC>CCC', 'CTC>CGC',
                   'CTG>CAG', 'CTG>CCG', 'CTG>CGG', 'CTT>CAT', 'CTT>CCT', 'CTT>CGT', 'GAA>GCA', 'GAA>GGA', 'GAA>GTA', 'GAC>GCC',
                   'GAC>GGC', 'GAC>GTC', 'GAG>GCG', 'GAG>GGG', 'GAG>GTG', 'GAT>GCT', 'GAT>GGT', 'GAT>GTT', 'GCA>GAA', 'GCA>GGA',
                   'GCA>GTA', 'GCC>GAC', 'GCC>GGC', 'GCC>GTC', 'GCG>GAG', 'GCG>GGG', 'GCG>GTG', 'GCT>GAT', 'GCT>GGT', 'GCT>GTT',
                   'GGA>GAA', 'GGA>GCA', 'GGA>GTA', 'GGC>GAC', 'GGC>GCC', 'GGC>GTC', 'GGG>GAG', 'GGG>GCG', 'GGG>GTG', 'GGT>GAT',
                   'GGT>GCT', 'GGT>GTT', 'GTA>GAA', 'GTA>GCA', 'GTA>GGA', 'GTC>GAC', 'GTC>GCC', 'GTC>GGC', 'GTG>GAG', 'GTG>GCG',
                   'GTG>GGG', 'GTT>GAT', 'GTT>GCT', 'GTT>GGT', 'TAA>TCA', 'TAA>TGA', 'TAA>TTA', 'TAC>TCC', 'TAC>TGC', 'TAC>TTC',
                   'TAG>TCG', 'TAG>TGG', 'TAG>TTG', 'TAT>TCT', 'TAT>TGT', 'TAT>TTT', 'TCA>TAA', 'TCA>TGA', 'TCA>TTA', 'TCC>TAC',
                   'TCC>TGC', 'TCC>TTC', 'TCG>TAG', 'TCG>TGG', 'TCG>TTG', 'TCT>TAT', 'TCT>TGT', 'TCT>TTT', 'TGA>TAA', 'TGA>TCA',
                   'TGA>TTA', 'TGC>TAC', 'TGC>TCC', 'TGC>TTC', 'TGG>TAG', 'TGG>TCG', 'TGG>TTG', 'TGT>TAT', 'TGT>TCT', 'TGT>TTT',
                   'TTA>TAA', 'TTA>TCA', 'TTA>TGA', 'TTC>TAC', 'TTC>TCC', 'TTC>TGC', 'TTG>TAG', 'TTG>TCG', 'TTG>TGG', 'TTT>TAT',
                   'TTT>TCT', 'TTT>TGT']
    dummy = pd.DataFrame({'from_to_context':dummy_list})

    df['counter'] = df.index
    data = df[['from_to_context', 'counter']].groupby(['from_to_context'], as_index=False).count()

    data = dummy.merge(data, how='left', on='from_to_context')
    data['counter'] = data['counter'].fillna(0)

    data['order1'] = data['from_to_context'].astype(str).str[1] + '>' + data['from_to_context'].astype(str).str[5]
    sorter1 = ['C>T', 'G>A', 'G>T', 'G>C', 'C>A', 'T>A', 'T>C', 'A>G', 'T>G', 'C>G', 'A>C', 'A>T']
    data.order1 = data.order1.astype("category")
    data.order1 = data.order1.cat.set_categories(sorter1)

    data['order2'] = data['from_to_context'].astype(str).str[0]
    sorter2 = ['C', 'T', 'G', 'A']
    data.order2 = data.order2.astype("category")
    data.order2 = data.order2.cat.set_categories(sorter2)

    data['order3'] = data['from_to_context'].astype(str).str[2]
    sorter3 = ['C', 'T', 'G', 'A']
    data.order3 = data.order3.astype("category")
    data.order3 = data.order3.cat.set_categories(sorter3)

    data = data.sort_values(["order1", "order2", "order3"])
    if replace_t_u:
        data['order1'] = data['order1'].str.replace('T', 'U')

    n_mut = int(sum(data['counter']))
    max_mut = int(max(data['counter']))
    data['counter'] = data['counter'] / n_mut
    if replace_t_u:
        data.loc[data['order1'].isin(['U>C', 'A>G', 'U>G', 'C>G', 'A>C', 'A>U']), 'counter'] = data['counter'] * -1
        color_mapping12 = {
            "U>G": "deepskyblue",
            "G>U": "deepskyblue",
            "U>A": "black",
            "A>U": "black",
            "C>U": "red",
            "U>C": "red",
            "C>A": "silver",
            "A>C": "silver",
            "G>A": "yellowgreen",
            "A>G": "yellowgreen",
            "G>C": "pink",
            "C>G": "pink",
        }
    else:
        data.loc[data['order1'].isin(['T>C', 'A>G', 'T>G', 'C>G', 'A>C', 'A>T']), 'counter'] = data['counter'] * -1
        color_mapping12 = {
            "T>G": "deepskyblue",
            "G>T": "deepskyblue",
            "T>A": "black",
            "A>T": "black",
            "C>T": "red",
            "T>C": "red",
            "C>A": "silver",
            "A>C": "silver",
            "G>A": "yellowgreen",
            "A>G": "yellowgreen",
            "G>C": "pink",
            "C>G": "pink",
        }


    data['color'] = data['order1']
    data = data.replace({"color": color_mapping12})
    if replace_t_u:
        data1 = data[data['order1'].isin(['C>U', 'G>A', 'G>U', 'G>C', 'C>A', 'U>A'])]
        data2 = data[data['order1'].isin(['U>C', 'A>G', 'U>G', 'C>G', 'A>C', 'A>U'])]
    else:
        data1 = data[data['order1'].isin(['C>T', 'G>A', 'G>T', 'G>C', 'C>A', 'T>A'])]
        data2 = data[data['order1'].isin(['T>C', 'A>G', 'T>G', 'C>G', 'A>C', 'A>T'])]

    fig, axes = plt.subplots(figsize=(20, 10), nrows=2)
    fig.subplots_adjust(hspace=0)

    axes[0].bar(data1['from_to_context'], data1['counter'], align='center', color=data1.color)

    axes[1].bar(data2['from_to_context'], data2['counter'], align='center', color=data2.color)

    axes[0].get_xaxis().set_ticks([])
    axes[1].get_xaxis().set_ticks([])

    if replace_t_u:
        text1 = ['C>U', 'G>A', 'G>U', 'G>C', 'C>A', 'U>A']
        text2 = ['U>C', 'A>G', 'U>G', 'C>G', 'A>C', 'A>U']
    else:
        text1 = ['C>T', 'G>A', 'G>T', 'G>C', 'C>A', 'T>A']
        text2 = ['T>C', 'A>G', 'T>G', 'C>G', 'A>C', 'A>T']
    color_text = ["red", 'yellowgreen', 'deepskyblue', 'pink', 'silver', 'black']

    pol = 8
    for text_num in range(len(text1)):
        plt.text(pol, max_mut / n_mut + 0.005, text1[text_num], fontsize=30, color=color_text[text_num])
        pol += 16

    pol = 8
    for text_num in range(len(text2)):
        plt.text(pol, -1 * (max_mut / n_mut) - 0.01, text2[text_num], fontsize=30, color=color_text[text_num])
        pol += 16

    axes[0].yaxis.set_tick_params(labelsize=20)
    axes[1].yaxis.set_tick_params(labelsize=20)

    if add_n:
        fig.suptitle('{} (n={})'.format(title, n_mut),
                     fontsize=30, horizontalalignment='left', x=0.15, y=1.1)
    else:
        fig.suptitle(title,
                     fontsize=30, horizontalalignment='left', x=0.15, y=1.1)

    plt.savefig(save_path, bbox_inches = 'tight')
