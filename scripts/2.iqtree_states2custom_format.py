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

# path_to_states = "../data/interim/iqtree_runs/brun3/anc_kg.state"
# path_to_leaves = "../data/interim/leaves_birds_states.tsv"
# out = "../data/interim/anc_kg_states_birds.tsv"

adtype = {
    "Site": np.int32,
    "p_A":  np.float32, "p_C":  np.float32,
    "p_G":  np.float32, "p_T":  np.float32,
}


def read_one_leaf(path_to_leaves):
    leaves = pd.read_csv(path_to_leaves, sep="\t", dtype=adtype)
    one_leaf = leaves[leaves.Node == leaves.Node.sample().values[0]].sort_values([
        "Part", "Site"])
    return one_leaf


def read_anc_states(path_to_states):
    anc = pd.read_csv(path_to_states, sep="\t", comment='#', dtype=adtype)
    anc = anc.sort_values(["Node", "Site"])
    return anc


@click.command("anc formatter", help="reformat anc states from iqtree output to custom states table")
@click.option("--anc", "path_to_states", required=True, type=click.Path(True), help="path to anc states from iqtree output")
@click.option("--leaves", "path_to_leaves", required=True, type=click.Path(True), help="path to leaves table of custom format")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output reformatted states file (tsv)")
def main(path_to_states, path_to_leaves, out):
    anc = read_anc_states(path_to_states)
    print("Ancestors read", file=sys.stderr)
    one_leaf = read_one_leaf(path_to_leaves)
    print("One leaf read", file=sys.stderr)
    replication_factor = anc.shape[0] / one_leaf.shape[0]
    assert replication_factor == int(replication_factor)
    replication_factor = int(replication_factor)

    anc["Part"] = np.tile(one_leaf.Part.values, replication_factor)
    anc["Site"] = np.tile(one_leaf.Site.values, replication_factor)
    anc = anc[["Node", "Part", "Site", "State", "p_A", "p_C", "p_G", "p_T"]]
    print("Modifications done\nWriting...", file=sys.stderr)
    anc.to_csv(out, sep="\t", index=None)


if __name__ == "__main__":
    main()
