#!/usr/bin/env python3

import os

import click
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from mutspec_utils.annotation import calculate_mutspec, rev_comp, transcriptor
from mutspec_utils.constants import possible_sbs192

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
colors12 = [color_mapping12[sbs] for sbs in sbs12_ordered]


kk_lbls = "A>C A>G A>T C>T G>C G>T".split()
cosmic_lbls = "C>A C>G C>T T>A T>C T>G".split()

df = pd.DataFrame({"sbs": possible_sbs192})
df["sbs_base"] = df["sbs"].str.slice(2, 5)
df["sbs_base_revcomp"] = df["sbs_base"].str.translate(transcriptor)
df["sbs_revcomp"] = df["sbs"].apply(rev_comp)
df["is_cosmic"] = df["sbs_base"].isin(cosmic_lbls)
df["is_kk"] = df["sbs_base"].isin(kk_lbls)
df["sbs_base_for_sorting_kp"] = df.apply(
    lambda x: x.sbs_base + "1" if x.is_cosmic else x.sbs_base_revcomp + "2", axis=1)
df["sbs_for_ordering_kk"] = df.apply(lambda x: x.sbs if x.is_kk else x.sbs_revcomp, axis=1)
df["sbs_for_ordering_kp"] = df.apply(lambda x: x.sbs if x.is_cosmic else x.sbs_revcomp, axis=1)

ordered_sbs192_kp = list(df.sort_values(["sbs_base_for_sorting_kp", "sbs_for_ordering_kp"]).sbs.values)
ordered_sbs192_kk = list(df.sort_values(["sbs_base", "sbs_for_ordering_kk"]).sbs.values)
del df


def __prepare_nice_labels(ordered_sbs192):
    _nice_order = []
    prev = None
    for sbs in ordered_sbs192:
        if prev is not None and sbs[2:5] != prev[2:5]:
            _nice_order.append("")
        # _nice_order.append(sbs[2] + sbs[4] + ": " + sbs[0] + sbs[2] + sbs[-1])
        _nice_order.append(sbs)
        prev = sbs
    return _nice_order


def plot_mutspec12bar(mutspec: pd.DataFrame, ylabel="MutSpec", title="Full mutational spectrum", show=True, savepath=None):
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax = sns.barplot(x="Mut", y=ylabel, data=mutspec, order=sbs12_ordered, ax=fig.gca())

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


def plot_mutspec12vio(mutspec: pd.DataFrame, ylabel="MutSpec", title="Full mutational spectrum", show=True, savepath=None):
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax = sns.violinplot(x="Mut", y=ylabel, data=mutspec, order=sbs12_ordered, ax=fig.gca())

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


def plot_mutspec192(mutspec192: pd.DataFrame, ylabel="MutSpec", title="Mutational spectrum", show=True, figsize=(24, 10), filepath=None):
    """
    Plot barblot of given mutational spectrum calculated from single nucleotide substitutions

    Arguments
    ---------
    mutspec192: pd.DataFrame
        table, containing 192 component mutational spectrum for one or many species, all substitutions must be presented in the table
    title: str, default = 'Mutational spectrum'
        Title on the plot
    filepath: str, default = None
        Path to output plot file. If None no images will be written
    """
    # TODO add checks of mutspec192
    ms192 = mutspec192.copy()
    ms192["MutBase"] = ms192.Mut.str.slice(2, 5)
    ms192["Context"] = ms192.Mut.str.get(0) + ms192.Mut.str.get(2) + ms192.Mut.str.get(-1)
    ms192["long_lbl"] = ms192.Mut.str.get(2) + ms192.Mut.str.get(4) + ": " + \
        ms192.Mut.str.get(0) + ms192.Mut.str.get(2) + ms192.Mut.str.get(-1)
    order = __prepare_nice_labels(ordered_sbs192_kp)

    # df = ms192.groupby(["MutBase", "Context"]).mean()
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.grid(axis="y", alpha=.7, linewidth=0.5)
    sns.barplot(
        x="Mut", y=ylabel, data=ms192,
        order=order, errwidth=1, ax=fig.gca(),
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    # map colors to bars
    width = 0.4
    shift = None
    for bar, sbs in zip(ax.patches, order):
        if not shift:
            shift = (bar.get_width() - width) / 2
        if len(sbs):
            bar.set_color(color_mapping12[sbs[2:5]])
            bar.set_alpha(alpha=0.9)
        bar.set_width(width)
        bar.set_x(bar.get_x() + shift)

    plt.xticks(rotation=90, fontsize=6)
    # labels = ['' for _ in ax.get_xticklabels()]
    # ax.set_xticklabels(labels)
    # __label_group_bar_table(ax, df)
    # fig.subplots_adjust(bottom=0.1 * df.index.nlevels)
    if filepath is not None:
        plt.savefig(filepath)
    if show:
        plt.show()
    else:
        plt.close()


def dump_expected(exp, path):
    exp_melted = exp.reset_index()
    exp_melted["Label"] = exp_melted["Label"].where(exp_melted["Label"] != "ff", "syn4f")
    exp_melted.melt("Label", exp_melted.columns.values[1:], var_name="Mut", value_name="Count")\
        .sort_values(["Mut", "Label"])\
            .to_csv(path, sep="\t", index=None)


@click.command("MutSpec calculator", help="Calculate and visualize mutational spectra")
@click.option("-b", "--observed", "path_to_obs", type=click.Path(True),  show_default=True, help="Path to observed mutations table")
@click.option("-e", "--expected", "path_to_exp", type=click.Path(True), help="Path to expected mutations table")
@click.option("-o", '--outdir', type=click.Path(True), default=".", show_default=True, help="Path to output directory for files (must exist)")
@click.option("-l", '--label', default="", show_default=True, help="Label for files naming")
@click.option("-p", '--proba', "use_proba", is_flag=True, help="Use probabilities of mutations")
@click.option('--proba_min', default=0.1, show_default=True, help="Minimal mutation probability to consider in mutspec calculation. Used only with --use_proba")
@click.option('--exclude', default="OUTGRP,ROOT", show_default=True, help="Name of source nodes to exclude from mutations set. Use comma to pass several names")
@click.option('--syn', is_flag=True, help="")
@click.option('--syn4f', is_flag=True, help="")
@click.option('--mnum192', "mut_num_for_192", default=16, show_default=True, help="Number of mutation types (maximum 192) required to calculate and plot 192-component mutational spectra")
@click.option('--plot', is_flag=True, help="Plot spectra plots")
@click.option("-x", '--ext', "image_extension", default="pdf", show_default=True, type=click.Choice(['pdf', 'png', 'jpg'], case_sensitive=False), help="Images format to save")
def main(path_to_obs, path_to_exp, outdir, label, use_proba, proba_min, exclude, syn, syn4f, mut_num_for_192, plot, image_extension):
    if mut_num_for_192 > 192:
        raise RuntimeError("Number of mutation types must be less then 192, but passed {}".format(mut_num_for_192))

    label_codes, labels = [0], ["all"]
    if syn:
        label_codes.append(1)
        labels.append("syn")
    if syn4f:
        label_codes.append(2)
        labels.append("ff")

    exclude = exclude.split(",")

    path_to_united_exp = os.path.join(outdir, "mean_expexted_mutations_{}.tsv".format(label))
    path_to_ms12 = os.path.join(outdir, "ms12{}_{}.tsv")
    path_to_ms12plot = os.path.join(outdir, "ms12{}_{}.{}")
    path_to_ms192plot = os.path.join(outdir, "ms192{}_{}.{}")
    path_to_ms192 = os.path.join(outdir, "ms192{}_{}.tsv")

    obs = pd.read_csv(path_to_obs, sep="\t")
    exp_raw = pd.read_csv(path_to_exp, sep="\t")

    # exclude ROOT node because RAxML don't change input tree that contains one node more than it's need
    obs = obs[~obs.AltNode.isin(exclude)]
    
    plot_mutspec12_func = plot_mutspec12vio
    if "Replica" not in obs.columns:
        obs["Replica"] = 1
        plot_mutspec12_func = plot_mutspec12bar

    if use_proba:
        obs = obs[(obs.ProbaFull > proba_min)]

    exp = exp_raw.drop_duplicates().drop(["Node", "Gene"], axis=1).groupby("Label").mean()
    dump_expected(exp, path_to_united_exp)

    for lbl_code, lbl in zip(label_codes, labels):
        cur_obs_lbl = obs[obs.Label >= lbl_code]
        if not cur_obs_lbl.shape[0]:
            continue

        if lbl not in exp.index:
            continue

        cur_exp = exp.loc[lbl].to_dict()

        ms12_collection = []
        ms192_collection = []
        for replica in cur_obs_lbl.Replica.unique():
            cur_obs_lbl_repl = cur_obs_lbl[cur_obs_lbl["Replica"] == replica]
            if not cur_obs_lbl_repl.shape[0]:
                continue
            ms12 = calculate_mutspec(cur_obs_lbl_repl, cur_exp, use_context=False, use_proba=use_proba)
            ms12["Mut"] = ms12["Mut"].str.translate(transcriptor)
            ms12.drop("RawMutSpec", axis=1, inplace=True)
            ms12_collection.append(ms12)

            if cur_obs_lbl_repl.Mut.nunique() >= mut_num_for_192:
                ms192 = calculate_mutspec(cur_obs_lbl_repl, cur_exp, use_context=True, use_proba=use_proba)
                ms192["Mut"] = ms192["Mut"].apply(rev_comp)
                ms192.drop("RawMutSpec", axis=1, inplace=True)
        
        if ms12_collection:
            ms12 = pd.concat(ms12_collection)
            ms12.to_csv(path_to_ms12.format(lbl, label), sep="\t", index=None)
            if plot:
                plot_mutspec12_func(
                    ms12, title=f"{lbl} mutational spectrum", 
                    savepath=path_to_ms12plot.format(lbl, label, image_extension), 
                    show=False
                )

        if ms192_collection:
            ms192 = pd.concat(ms192_collection)
            ms192.to_csv(path_to_ms192.format(lbl, label), sep="\t", index=None)
            if plot:
                plot_mutspec192(
                    ms192, title=f"{lbl} mutational spectrum",
                    filepath=path_to_ms192plot.format(lbl, label, image_extension), 
                    show=False
                )


if __name__ == "__main__":
    main()
    # main("-b tmp/observed_mutations_iqtree.tsv -e tmp/expected_mutations.tsv -o tmp/ -l debug".split())
