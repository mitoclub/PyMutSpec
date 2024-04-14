#!/usr/bin/env python3

import os
import sys
from functools import partial

import click
import pandas as pd

from pymutspec.annotation import calculate_mutspec
from pymutspec.constants import possible_sbs12, possible_sbs192
from pymutspec.draw import plot_mutspec12, plot_mutspec192


def save_tsv(df: pd.DataFrame, path):
    df.to_csv(path, sep="\t", float_format='%g', index=False)


def dump_expected(exp, path):
    exp_melted = exp.reset_index()
    exp_melted["Label"] = exp_melted["Label"].where(exp_melted["Label"] != "ff", "syn4f")
    exp_melted = exp_melted.melt(
        "Label", exp_melted.columns.values[1:], var_name="Mut", value_name="Count"
    ).sort_values(["Mut", "Label"])
    save_tsv(exp_melted, path)


@click.command("MutSpec calculator", help="Calculate and visualize mutational spectra")
@click.option("-b", "--observed", "path_to_obs", type=click.Path(True),  show_default=True, help="Path to observed mutations table")
@click.option("-e", "--expected", "path_to_exp", type=click.Path(True), help="Path to expected mutations table")
@click.option("-o", '--outdir', type=click.Path(True), default=".", show_default=True, help="Path to output directory for files (must exist)")
@click.option("-l", '--label', default=None, help="Label for files naming. By default no label")
@click.option("-p", '--proba', "use_proba", is_flag=True, help="Use probabilities of mutations")
@click.option('--proba_min', default=0.3, show_default=True, help="LEGACY Minimal mutation probability to consider in spectra calculation. Used only with --use_proba")
@click.option('--proba_cutoff', default=None, type=float, show_default=True, help="Minimal mutation probability to consider in spectra calculation. Used only with --use_proba")
@click.option('--exclude', default="OUTGRP,ROOT", show_default=True, help="Name of source nodes to exclude from mutations set. Use comma to pass several names")
@click.option('--all', 'all_muts', is_flag=True, help="Calculate and plot spectra for all mutations; "
                                                      "default if not specified at least one of --all, --syn, --syn4f")
@click.option('--syn', is_flag=True, help="Calculate and plot spectra for synonymous mutations")
@click.option('--syn4f', is_flag=True, help="Calculate and plot spectra for synonymous mutations in fourfols positions")
@click.option('--nonsyn', is_flag=True, help="Calculate and plot spectra for non-synonymous mutations")
@click.option('--mnum192', type=click.IntRange(0, 192), default=16, show_default=True, help="Number of mutation types (maximum 192) required to calculate and plot 192-component mutational spectra")
@click.option('--substract12', "path_to_substract12",   type=click.Path(True), default=None, help="Mutational spectrum that will be substracted from calculated spectra 12 component")
@click.option('--substract192', "path_to_substract192", type=click.Path(True), default=None, help="Mutational spectrum that will be substracted from calculated spectra 192 component")
@click.option('--branches', is_flag=True, help="Calculate tree branch specific spectra")
@click.option('--subset', default="full", show_default=True, type=click.Choice(["full", "internal", "terminal"]), help="Used subset of tree branches")
@click.option('--plot', is_flag=True, help="Plot spectra plots")
@click.option("-x", '--ext', "image_extension", default="pdf", show_default=True, type=click.Choice(['pdf', 'png', 'jpg'], case_sensitive=False), help="Images format to save")
def main(
        path_to_obs, path_to_exp, outdir, label, 
        use_proba, proba_min, proba_cutoff, exclude, 
        all_muts, syn, syn4f, nonsyn, mnum192, 
        path_to_substract12, path_to_substract192, 
        branches, subset, plot, image_extension,
    ):
    proba_cutoff = proba_cutoff or proba_min

    if mnum192 > 192:
        raise RuntimeError("Number of mutation types must be less then 192, but passed {}".format(mnum192))

    lbl_ids, lbls = [], []
    if all_muts:
        lbl_ids.append(0)
        lbls.append("all")
    if syn:
        lbl_ids.append(1)
        lbls.append("syn")
    if syn4f:
        lbl_ids.append(2)
        lbls.append("ff")
    if nonsyn:
        lbl_ids.append(None) # ???
        lbls.append("nonsyn")
    if len(lbls) == 0:
        lbl_ids, lbls = [0], ["all"]

    exclude = exclude.split(",")
    
    if label is None:
        label = "" if subset == "full" else "_" + subset
    else:
        label = "_" + label if subset == "full" else "_{}_{}".format(subset, label)

    path_to_ms12      = os.path.join(outdir, "ms12{}{}.tsv")
    path_to_ms12plot  = os.path.join(outdir, "ms12{}{}.{}")
    path_to_ms192plot = os.path.join(outdir, "ms192{}{}.{}")
    path_to_ms192     = os.path.join(outdir, "ms192{}{}.tsv")

    substract12  = pd.read_csv(path_to_substract12,  sep="\t") if path_to_substract12  is not None else None
    substract192 = pd.read_csv(path_to_substract192, sep="\t") if path_to_substract192 is not None else None

    obs = pd.read_csv(path_to_obs, sep="\t")
    # exclude ROOT node because RAxML don't change input tree that contains one node more than it's need
    obs = obs[~obs.AltNode.isin(exclude)]
    if subset == "internal":
        # TODO replace .startswith("Node") to more common method: extract leaves names with PhyloTree and ...
        obs = obs[obs.AltNode.str.startswith("Node")]
    elif subset == "terminal":
        obs = obs[~obs.AltNode.str.startswith("Node")]

    exp_raw = pd.read_csv(path_to_exp, sep="\t")
    # if all substitutions are columns
    if len(set(possible_sbs12 + possible_sbs192).difference(exp_raw.columns)) == 0:
        exp_mean = exp_raw.drop_duplicates().drop(["Node", "Gene"], axis=1, errors="ignore").groupby("Label").mean()
        path_to_united_exp = os.path.join(outdir, "mean_expexted_mutations{}.tsv".format(label))
        dump_expected(exp_mean, path_to_united_exp)
        exp_freqs = exp_raw
        del exp_raw
    elif list(exp_raw.columns) == ["Label", "Mut", "Count"]:
        if branches:
            raise ValueError("For branch specific spectra expected mutations required for every internal tree node")

        exp_mean = exp_raw.pivot("Label", "Mut", "Count")
        exp_freqs = None
    else:
        raise RuntimeError("Expected another columns in the table {}".format(path_to_exp))
    
    plot_mutspec12_func = partial(plot_mutspec12, style="box")
    if "Replica" not in obs.columns:
        obs["Replica"] = 1
        plot_mutspec12_func = partial(plot_mutspec12, style="bar")
    if substract12 is not None:
        plot_mutspec12_func = partial(plot_mutspec12, style="bar")

    if use_proba:
        obs = obs[(obs.ProbaFull > proba_cutoff)]

    for lbl_id, lbl in zip(lbl_ids, lbls):
        if lbl == "nonsyn":
            cur_obs = obs[obs.Label == 0]
        else:
            cur_obs = obs[obs.Label >= lbl_id]

        if not cur_obs.shape[0] or lbl not in exp_mean.index:
            continue

        cur_exp = exp_mean.loc[lbl].to_dict()

        ms12_collection  = []
        ms192_collection = []
        for replica in cur_obs["Replica"].unique():
            cur_obs_repl = cur_obs[cur_obs["Replica"] == replica]
            if not cur_obs_repl.shape[0]:
                continue
            ms12 = calculate_mutspec(cur_obs_repl, cur_exp, use_context=False, use_proba=use_proba)
            ms12_collection.append(ms12)

            if cur_obs_repl.Mut.nunique() >= mnum192:
                ms192 = calculate_mutspec(cur_obs_repl, cur_exp, use_context=True, use_proba=use_proba)
                ms192_collection.append(ms192)
                
        
        if ms12_collection:
            ms12 = pd.concat(ms12_collection)
            save_tsv(ms12, path_to_ms12.format(lbl, label))
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
            save_tsv(ms192, path_to_ms192.format(lbl, label))
            if plot:
                if substract192 is None:
                    plot_mutspec192(
                        ms192[(ms192["ObsNum"] > 1) & (ms192["ExpNum"] > 1)], 
                        title=f"{lbl} mutational spectrum",
                        savepath=path_to_ms192plot.format(lbl, label, image_extension), 
                        show=False,
                    )
                else:
                    # TODO don't show useless bars
                    substract192.loc[(substract192["ObsNum"] <= 1) | (substract192["ExpNum"] <= 1), "MutSpec"] = 0.
                    substract192 = substract192.rename(columns={"MutSpec": "MutSpec_obs"})[["Mut", "MutSpec_obs"]]
                    ms192 = ms192.rename(columns={"MutSpec": "MutSpec_exp"}).merge(substract192, on="Mut")
                    ms192["MutSpec"] = ms192["MutSpec_obs"] - ms192["MutSpec_exp"]
                    plot_mutspec192(
                        ms192, 
                        title=f"{lbl} mutational spectrum difference\n(reconstructed - simulated (obs - exp))",
                        savepath=path_to_ms192plot.format(lbl, label, image_extension), 
                        show=False,
                    )
        if branches:
            # prepare branch-specific frequencies (sbs12 + sbs192)
            exp_freqs_lbl = exp_freqs[exp_freqs.Label == lbl].drop(["Gene", "Label"], axis=1).set_index('Node')
            assert exp_freqs_lbl.shape[1] == 192 + 12, f'prepared expected freqs absent some substitutions; freqs shape: {exp_freqs_lbl.shape}'

            if substract12 is not None or substract192 is not None:
                print("Ignoring substract agruments, substraction is not available for branch spectra", file=sys.stderr)
            
            branch_mutspec12, branch_mutspec192 = [], []
            for alt_node in cur_obs["AltNode"].unique():
                branch_obs = cur_obs[cur_obs["AltNode"] == alt_node]
                ref_node = branch_obs['RefNode'].iloc[0]
                branch_exp = exp_freqs_lbl.loc[ref_node].to_dict()
                if len(branch_obs) > 1:
                    # at least 2 mutations on branch
                    ms12 = calculate_mutspec(branch_obs, branch_exp, use_context=False, use_proba=use_proba)\
                        .assign(RefNode=ref_node, AltNode=alt_node, )
                    branch_mutspec12.append(ms12)

                if branch_obs.Mut.nunique() >= mnum192 and branch_obs.shape[0] > 10:
                    # at least 10 mutations on branch
                    ms192 = calculate_mutspec(branch_obs, branch_exp, use_context=True, use_proba=use_proba)\
                        .assign(RefNode=ref_node, AltNode=alt_node, )
                    branch_mutspec192.append(ms192)
            
            if len(branch_mutspec12) > 0:
                branch_mutspec12df  = pd.concat(branch_mutspec12)
                save_tsv(branch_mutspec12df, path_to_ms12.format(lbl, label).replace(".tsv", "_branches.tsv"))

            if len(branch_mutspec192) > 0:
                branch_mutspec192df = pd.concat(branch_mutspec192)
                save_tsv(branch_mutspec192df, path_to_ms192.format(lbl, label).replace(".tsv", "_branches.tsv"))


if __name__ == "__main__":
    main()
    # pt = '/home/kpotoh/nemu-pipeline/example_mam_run/outdir/mutspec_tables'
    # main(f"-b {pt}/observed_mutations_iqtree.tsv -e {pt}/expected_freqs.tsv -o {pt}/test_branches -l debug --branches --syn".split())
