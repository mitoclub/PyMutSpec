#!/usr/bin/env python3

import re
import sys

import click
from Bio import SeqIO


def pattern_selector(records, patterns: list):
    patterns = [re.compile(p) for p in patterns]
    for rec in records:
        for p in patterns:
            if p.search(rec.description):
                yield rec
                break

def read_patterns_from_file(fname: str) -> list:
    p = []
    with open(fname) as fin:
        for line in fin:
            p.append(line.strip())
    return p


@click.command()
@click.argument("aln", type=click.Path(True))
@click.argument("outfile", type=click.Path(writable=True), required=False, default=None)
@click.option("-p", "--pattern", multiple=True, default=None,
                help="Pattern used to select records by header; \
                    must be interpretable by re module; \
                        several patterns could be used")
@click.option("-f", "--file", type=click.Path(True), default=None, help="Path to file containing patterns (one pattern per line)")
@click.option("--fasta", "fasta_format", is_flag=True, )
@click.option("--phylip", "phylip_format", is_flag=True, )
def main(aln, outfile, pattern, file, fasta_format, phylip_format):
    """
    Select sequences from fasta file by header pattern

    if no OUTFILE passed will print sequences to stdout
    """
    if (len(pattern) and file is None) or (len(pattern) > 0 and file is not None):
        raise RuntimeError("Either File or Pattern must be selected as input")
    
    if fasta_format:
        fmt = "fasta"
        fmt_out = "fasta-2line"
    elif phylip_format:
        fmt = fmt_out = "phylip"
    else:
        raise RuntimeError("Choose the alignment format: --fasta or --phylip")
    
    if len(pattern) == 0:
        pattern = read_patterns_from_file(file)

    records = SeqIO.parse(aln, format=fmt)
    selected = pattern_selector(records, pattern)
    outfile = outfile or sys.stdout

    SeqIO.write(selected, outfile, fmt_out)


if __name__ == "__main__":
    main()
