#!/usr/bin/env python3

import re
import sys

import click
from Bio import SeqIO


def pattern_selector(records, patterns: list, invert_match=False):
    patterns = [re.compile(p) for p in patterns]
    is_match = False
    for rec in records:
        for p in patterns:
            if p.search(rec.description):
                is_match = True
                break
        if not invert_match and is_match:
            yield rec
        elif invert_match and not is_match:
            yield rec

        is_match = False


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
@click.option("-v", "--invert-match", is_flag=True, help="Invert the sense of matching, to select non-matching lines")
@click.option("--fmt", type=click.Choice(["fasta", "phylip"]),  default="fasta", show_default=True, help="Format of input alignment")
@click.option("--fasta", "fasta_format", is_flag=True, help="Legacy, use --fmt option")
@click.option("--phylip", "phylip_format", is_flag=True, help="Legacy, use --fmt option")
def main(aln, outfile, pattern, file, invert_match, fmt, fasta_format, phylip_format):
    """
    Select sequences from fasta file by header pattern

    if no OUTFILE passed will print sequences to stdout
    """
    if (len(pattern) == 0 and file is None) or (len(pattern) > 0 and file is not None):
        raise RuntimeError(f"Either File or Pattern must be selected as input, passed pattern ({pattern}) and file ({file})")
    
    if phylip_format or fmt == "phylip":
        fmt = fmt_out = "phylip"
    elif fasta_format or fmt == "fasta":
        fmt = "fasta"
        fmt_out = "fasta-2line"
    else:
        raise RuntimeError()
    
    if len(pattern) == 0:
        pattern = read_patterns_from_file(file)

    records = SeqIO.parse(aln, format=fmt)
    selected = pattern_selector(records, pattern, invert_match)
    outfile = outfile or sys.stdout

    SeqIO.write(selected, outfile, fmt_out)


if __name__ == "__main__":
    main()
