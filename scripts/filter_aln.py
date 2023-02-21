#!/usr/bin/env python3

import click
from Bio import SeqIO

FORMAT = "phylip"


@click.command("Alignment Filtrator", help="Filter out sequences from alignment according to the tree nodes list")
@click.option("-a", "--aln", type=click.Path(True), help="Path to input alignment that will be filtered")
@click.option("-l", "--leaves", type=click.Path(True),  show_default=True, help="Path to tree leaves list, one node name per one line (txt)")
@click.option("-o", "--out_aln", type=click.Path(writable=True), help="Path to output filtered alignment")
def main(leaves, aln, out_aln):
    seqs = SeqIO.parse(aln, FORMAT)
    nodes = set()
    with open(leaves) as fin:
        for line in fin:
            nodes.add(line.strip())

    out_seqs = []
    for r in seqs:
        if r.id in nodes:
            out_seqs.append(r)

    SeqIO.write(out_seqs, out_aln, FORMAT)


if __name__ == "__main__":
    main()
