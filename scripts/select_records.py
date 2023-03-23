#!/usr/bin/env python3

import re
import sys

import click
from Bio import SeqIO


def pattern_selector(fasta, patterns: list):
    patterns = [re.compile(p) for p in patterns]
    for rec in SeqIO.parse(fasta, "fasta"):
        for p in patterns:
            if p.search(rec.description):
                yield rec
                break


@click.command()
@click.argument("fasta", type=click.Path(True))
@click.argument("outfile", type=click.Path(writable=True), required=False, default=None)
@click.option("-p", "--pattern", multiple=True, required=True, 
                help="Pattern used to select records by header; \
                    must be interpretable by re module; \
                        several patterns could be used")
def main(fasta, outfile, pattern):
    """
    Select sequences from fasta file by header pattern

    if no OUTFILE passed will print sequences to stdout
    """
    selected = pattern_selector(fasta, pattern)
    outfile = outfile or sys.stdout
    SeqIO.write(selected, outfile, "fasta-2line")


if __name__ == "__main__":
    main()
