import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mutspec_utils.draw import plot_mutspec192kk, plot_mutspec192, plot_mutspec12
from mutspec_utils.constants import possible_sbs12, possible_sbs192
from mutspec_utils.annotation import calculate_mutspec

OUTGRP = "OUTGRP"
LABELS = [0, 1, 2]
MUT_NUM_FOR_192 = 150

# plot_mutspec192kk(df, )
# plot_mutspec192(df)

# @click.command("MutSpec calculator", help="TODO")
# @click.option("--db_path", "path_to_db", type=click.Path(writable=True), default="/tmp/states.db", show_default=True, help="Path to database with states. Use only with --write_db")
# @click.option("--rewrite_db",  is_flag=True, default=False, help="Rewrite existing states database. Use only with --write_db")
# @click.option('-f', '--force', is_flag=True, help="Rewrite existing output directory")
# @click.option('-v', '--verbose', "verbosity", count=True, help="Verbosity level = DEBUG")
def main():
    path_to_obs = "./tmp/mutations.tsv"
    path_to_exp = "./tmp/expected_mutations.tsv"
    path_to_ms12 = "./tmp/ms12{}.tsv"
    path_to_ms12plot = "./tmp/ms12{}.png"
    path_to_ms192plot = "./tmp/ms192{}.png"
    path_to_ms192 = "./tmp/ms192{}.tsv"

    obs = pd.read_csv(path_to_obs, sep="\t")
    exp_raw = pd.read_csv(path_to_exp, sep="\t")

    obs = obs[obs.AltNode != OUTGRP]
    exp = exp_raw.drop_duplicates().drop(["Node", "Gene"], axis=1).groupby("Label").mean()

    for lbl_code in LABELS:
        if lbl_code == 0:
            lbl = "all"
        elif lbl_code == 1:
            lbl = "syn"
        elif lbl_code == 2:
            lbl = "ff"

        cur_obs = obs[obs.Label >= lbl_code]
        cur_exp = exp.loc[lbl].to_dict()
        
        ms12 = calculate_mutspec(cur_obs, cur_exp, use_context=False, use_proba=True)
        ms12.to_csv(path_to_ms12.format(lbl), sep="\t", index=None)
        plot_mutspec12(ms12, savepath=path_to_ms12plot.format(lbl), show=False)
        if cur_obs.shape[0] > MUT_NUM_FOR_192:
            ms192 = calculate_mutspec(cur_obs, cur_exp, use_context=True, use_proba=True)
            ms192.to_csv(path_to_ms192.format(lbl), sep="\t", index=None)
            plot_mutspec192(ms192, filepath=path_to_ms192plot.format(lbl))


if __name__ == "__main__":
    main()
