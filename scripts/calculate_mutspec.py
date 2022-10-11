import os

import click
import pandas as pd

from mutspec_utils.draw import plot_mutspec192, plot_mutspec12
from mutspec_utils.annotation import calculate_mutspec

OUTGRP = "OUTGRP"
MUT_NUM_FOR_192 = 150
LABELS = [0, 1, 2]


@click.command("MutSpec calculator", help="Calculate and visualize mutational spectra")
@click.option("-b", "--observed", "path_to_obs", type=click.Path(True),  show_default=True, help="Path to observed mutations table")
@click.option("-e", "--expected", "path_to_exp", type=click.Path(True), help="Path to expected mutations table")
@click.option("-o", '--outdir', type=click.Path(True), default=".", show_default=True, help="Path to output directory for files (must exist)")
@click.option('--outgrp', default=OUTGRP, show_default=True, help="Name of outgroup node in the tree to exclude from mutations set")
@click.option('--mnum192', default=MUT_NUM_FOR_192, show_default=True, help="Number of mutations required to calculate and plot 192-component mutational spectra")
def main(path_to_obs, path_to_exp, outdir, outgrp, mut_num_for_192):
    path_to_united_exp = os.path.join(outdir, "mean_expexted_mutations.tsv")
    path_to_ms12 = os.path.join(outdir, "ms12{}.tsv")
    path_to_ms12plot = os.path.join(outdir, "ms12{}.png")
    path_to_ms192plot = os.path.join(outdir, "ms192{}.png")
    path_to_ms192 = os.path.join(outdir, "ms192{}.tsv")

    obs = pd.read_csv(path_to_obs, sep="\t")
    exp_raw = pd.read_csv(path_to_exp, sep="\t")

    obs = obs[obs.AltNode != outgrp]
    exp = exp_raw.drop_duplicates().drop(["Node", "Gene"], axis=1).groupby("Label").mean()
    exp.reset_index().to_csv(path_to_united_exp, sep="\t", index=None)

    for lbl_code in LABELS:
        if lbl_code == 0:
            lbl = "all"
        elif lbl_code == 1:
            lbl = "syn"
        elif lbl_code == 2:
            lbl = "ff"

        cur_obs = obs[obs.Label >= lbl_code]
        if cur_obs.shape[0]:
            cur_exp = exp.loc[lbl].to_dict()

            ms12 = calculate_mutspec(cur_obs, cur_exp, use_context=False, use_proba=True)
            ms12.to_csv(path_to_ms12.format(lbl), sep="\t", index=None)
            plot_mutspec12(ms12, savepath=path_to_ms12plot.format(lbl), show=False)
            if cur_obs.shape[0] > mut_num_for_192:
                ms192 = calculate_mutspec(cur_obs, cur_exp, use_context=True, use_proba=True)
                ms192.to_csv(path_to_ms192.format(lbl), sep="\t", index=None)
                plot_mutspec192(ms192, filepath=path_to_ms192plot.format(lbl))


if __name__ == "__main__":
    main()
