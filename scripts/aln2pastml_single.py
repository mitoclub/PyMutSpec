#!/usr/bin/env python3

import click
import pandas as pd
from Bio import SeqIO

GAP_CUTOFF = 0.1


def aln2pastml(filepath):
    data = []
    fasta = SeqIO.parse(filepath, "fasta")
    for rec in fasta:
        node = rec.name
        seq = str(rec.seq)
        for site, state in enumerate(seq, 1):
            pos_data = [node, site, state]
            data.append(pos_data)

    df = pd.DataFrame(data, columns="Node Site State".split()).sort_values(["Node", "Site"])
    pidf = df.pivot("Node", "Site", "State").reset_index()
    condition = (pidf == "-").sum(axis=0) < pidf.shape[0] * GAP_CUTOFF
    pidf = pidf.loc[:, condition]
    return pidf


@click.command("formatter", help="reformat alignment to states table")
@click.option("--aln", required=True, type=click.Path(True), help="path to file with single gene alignment")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(aln, out):
    df = aln2pastml(aln)
    df.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
    # main("data/example_nematoda/alignments_nematoda_clean", "data/example_nematoda/scheme_devilworm.nex", "data/example_nematoda/leaves")
