#!/usr/bin/env python3

import os
from functools import partial

import click
import pandas as pd

from pymutspec.annotation import calculate_mutspec, rev_comp, transcriptor
from pymutspec.draw import plot_mutspec12, plot_mutspec192


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
@click.option("-l", '--label', default=None, show_default=True, help="Label for files naming. By default no label")
@click.option("-p", '--proba', "use_proba", is_flag=True, help="Use probabilities of mutations")
@click.option('--proba_min', default=0.1, show_default=True, help="Minimal mutation probability to consider in spectra calculation. Used only with --use_proba")
@click.option('--exclude', default="OUTGRP,ROOT", show_default=True, help="Name of source nodes to exclude from mutations set. Use comma to pass several names")
@click.option('--all', 'all_muts', is_flag=True, help="Calculate and plot spectra for all mutations; "
                                                      "default if not specified at least one of --all, --syn, --syn4f")
@click.option('--syn', is_flag=True, help="Calculate and plot spectra for synonymous mutations")
@click.option('--syn4f', is_flag=True, help="Calculate and plot spectra for synonymous mutations in fourfols positions")
@click.option('--mnum192', "mut_num_for_192", default=16, show_default=True, help="Number of mutation types (maximum 192) required to calculate and plot 192-component mutational spectra")
@click.option('--substract12', "path_to_substract12",   type=click.Path(True), default=None, help="Mutational spectrum that will be substracted from calculated spectra 12 component")
@click.option('--substract192', "path_to_substract192", type=click.Path(True), default=None, help="Mutational spectrum that will be substracted from calculated spectra 192 component")
@click.option('--plot', is_flag=True, help="Plot spectra plots")
@click.option("-x", '--ext', "image_extension", default="pdf", show_default=True, type=click.Choice(['pdf', 'png', 'jpg'], case_sensitive=False), help="Images format to save")
def main(
    path_to_obs, path_to_exp, outdir, label, use_proba, proba_min, exclude, 
    all_muts, syn, syn4f, mut_num_for_192, path_to_substract12, path_to_substract192, 
    plot, image_extension,
    ):
    if mut_num_for_192 > 192:
        raise RuntimeError("Number of mutation types must be less then 192, but passed {}".format(mut_num_for_192))

    label_codes, labels = [], []
    if all_muts:
        label_codes.append(0)
        labels.append("all")
    if syn:
        label_codes.append(1)
        labels.append("syn")
    if syn4f:
        label_codes.append(2)
        labels.append("ff")
    if len(labels) == 0:
        label_codes, labels = [0], ["all"]

    exclude = exclude.split(",")
    
    label = "" if label is None else "_" + label
    path_to_ms12 = os.path.join(outdir, "ms12{}{}.tsv")
    path_to_ms12plot = os.path.join(outdir, "ms12{}{}.{}")
    path_to_ms192plot = os.path.join(outdir, "ms192{}{}.{}")
    path_to_ms192 = os.path.join(outdir, "ms192{}{}.tsv")

    substract12 = pd.read_csv(path_to_substract12, sep="\t") if path_to_substract12 is not None else None
    substract192 = pd.read_csv(path_to_substract192, sep="\t") if path_to_substract192 is not None else None

    obs = pd.read_csv(path_to_obs, sep="\t")
    # exclude ROOT node because RAxML don't change input tree that contains one node more than it's need
    obs = obs[~obs.AltNode.isin(exclude)]

    exp_raw = pd.read_csv(path_to_exp, sep="\t")
    # if random codons in columns
    if len(exp_raw.columns) > 192+12 and "A[T>C]A" in exp_raw.columns and "T[C>T]T" in exp_raw.columns:
        exp = exp_raw.drop_duplicates().drop(["Node", "Gene"], axis=1, errors="ignore").groupby("Label").mean()
        path_to_united_exp = os.path.join(outdir, "mean_expexted_mutations{}.tsv".format(label))
        dump_expected(exp, path_to_united_exp)
    elif list(exp_raw.columns) == ["Label", "Mut", "Count"]:
        exp = exp_raw.pivot("Label", "Mut", "Count")
    else:
        raise RuntimeError("Expected another columns in the table {}".format(path_to_exp))
    
    plot_mutspec12_func = partial(plot_mutspec12, style="box")
    if "Replica" not in obs.columns or substract12 is not None:
        obs["Replica"] = 1
        plot_mutspec12_func = partial(plot_mutspec12, style="bar")

    if use_proba:
        obs = obs[(obs.ProbaFull > proba_min)]

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
                ms192_collection.append(ms192)
                
        
        if ms12_collection:
            ms12 = pd.concat(ms12_collection)
            ms12.to_csv(path_to_ms12.format(lbl, label), sep="\t", index=None)
            if plot:
                if substract12 is None:
                    plot_mutspec12_func(
                        ms12, 
                        title=f"{lbl} mutational spectrum", 
                        savepath=path_to_ms12plot.format(lbl, label, image_extension), 
                        show=False,
                    )
                else:
                    ms12 = ms12.rename(columns={"MutSpec": "MutSpec_exp"})\
                        .merge(substract12.rename(columns={"MutSpec": "MutSpec_obs"})[["Mut", "MutSpec_obs"]], on="Mut")
                    ms12["MutSpec"] = ms12["MutSpec_obs"] - ms12["MutSpec_exp"]
                    plot_mutspec12_func(
                        ms12, 
                        title=f"{lbl} mutational spectrum difference\n(reconstructed - simulated (obs - exp))", 
                        savepath=path_to_ms12plot.format(lbl, label, image_extension), 
                        show=False,
                    )

        if ms192_collection:
            ms192 = pd.concat(ms192_collection)
            ms192.to_csv(path_to_ms192.format(lbl, label), sep="\t", index=None)
            if plot:
                if substract192 is None:
                    plot_mutspec192(
                        ms192, 
                        title=f"{lbl} mutational spectrum",
                        savepath=path_to_ms192plot.format(lbl, label, image_extension), 
                        show=False,
                        fontsize=7,
                    )
                else:
                    ms192 = ms192.rename(columns={"MutSpec": "MutSpec_exp"})\
                        .merge(substract192.rename(columns={"MutSpec": "MutSpec_obs"})[["Mut", "MutSpec_obs"]], on="Mut")
                    ms192["MutSpec"] = ms192["MutSpec_obs"] - ms192["MutSpec_exp"]
                    plot_mutspec192(
                        ms192, 
                        title=f"{lbl} mutational spectrum difference\n(reconstructed - simulated (obs - exp))",
                        savepath=path_to_ms192plot.format(lbl, label, image_extension), 
                        show=False,
                        fontsize=7,
                    )


if __name__ == "__main__":
    main()
    # main("-b tmp/observed_mutations_iqtree.tsv -e tmp/expected_mutations.tsv -o tmp/ -l debug".split())
