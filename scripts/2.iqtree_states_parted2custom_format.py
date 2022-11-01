#!/usr/bin/env python3

"""
Adapt anc table (part adding & site replacement)

- Part is gene index from scheme (but here used leaves annotation that include scheme)
- Site is gene site, not alignment
"""
import sys

import click
import numpy as np
import pandas as pd

from mutspec_utils.io import load_scheme

dtypes = {
    "Site": np.int32,
    "p_A":  np.float32, "p_C":  np.float32,
    "p_G":  np.float32, "p_T":  np.float32,
}


def read_anc_states(path_to_states):
    anc = pd.read_csv(path_to_states, sep="\t", comment='#', dtype=dtypes)
    anc = anc.sort_values(["Node", "Part", "Site"])
    return anc


@click.command("anc formatter", help="reformat anc states from iqtree output to custom states table")
@click.option("--anc", "path_to_states", required=True, type=click.Path(True), help="path to anc states from iqtree output")
@click.option("--scheme", "path_to_scheme", required=True, type=click.Path(True), help="path to scheme")
@click.option("--leaves", "path_to_leaves", required=False, type=click.Path(True), default=None, help="path to leaves table of custom format, used for QC")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output reformatted states file (tsv)")
def main(path_to_states, path_to_scheme, path_to_leaves, out):
    scheme_map = load_scheme(path_to_scheme)
    print("Scheme loaded", file=sys.stderr)

    anc = read_anc_states(path_to_states)
    print("Ancestors read", file=sys.stderr)

    anc["Part"] = anc["Part"].map(scheme_map).str.removesuffix(".fna")

    if path_to_leaves is not None:
        leaves = pd.read_csv(path_to_leaves, sep="\t", dtype=dtypes)
        dl = (leaves.groupby(["Part"]).Site.count() / leaves.Node.nunique()).to_dict()
        da = (anc.groupby(["Part"]).Site.count() / anc.Node.nunique()).to_dict()
        if dl != da:
            raise RuntimeError("Leaves states don't correspond to anc states")

    anc = anc[["Node", "Part", "Site", "State", "p_A", "p_C", "p_G", "p_T"]]
    print("Modifications done\nWriting...", file=sys.stderr)
    anc.to_csv(out, sep="\t", index=None)


if __name__ == "__main__":
    main()
    # main(
    #     "./data/example_nematoda/nematoda_anc_HKY_part/anc_HKY_part.state", 
    #     "./data/example_nematoda/scheme_devilworm.nex", 
    #     "./data/example_nematoda/leaves_states_nematoda.tsv", 
    #     "/tmp/trable.tsv"
    # )
