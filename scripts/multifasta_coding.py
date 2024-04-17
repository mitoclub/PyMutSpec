#!/usr/bin/env python3

import click
from Bio import SeqIO


def is_coded(aln):
    ncoded = noutgrp = 0
    for rec in aln:
        if rec.description == "OUTGRP":
            noutgrp += 1
        elif rec.description.startswith("RN_"):
            ncoded += 1
    return ncoded == len(aln) - 1 and noutgrp == 1


def is_outgroup(header: str, outgroup: str):
    is_equal = header == outgroup
    is_substring = outgroup in header
    return is_equal or is_substring


@click.command("Fasta header coder", help="")
@click.option("-a", "--alignment", type=click.Path(True), help="Fasta nucleotide alignment")
@click.option("-g", "--outgroup", help="Outrgoup header. Can be substring of header")
@click.option("-o", "--outfile", type=click.Path(writable=True), help="Output fasta file with coded headers")
@click.option("-m", "--outmap", type=click.Path(writable=True), help="Mapping of headers to coders (txt)")
def main(alignment, outgroup, outfile, outmap):
    aln = [r for r in SeqIO.parse(alignment, "fasta")]

    if is_coded(aln):
        print("Coding is not required, headers already coded")
    else:
        sp_map = dict()
        for i, rec in enumerate(aln, 1):
            header = rec.description
            coded_header = "OUTGRP" if is_outgroup(header, outgroup) else "RN_{}".format(i)
            sp_map[coded_header] = header

            rec.id = coded_header
            rec.name = coded_header
            rec.description = coded_header

        with open(outmap, "w") as fout:
            fout.write("\n".join(["{}\t{}".format(k, v) for k,v in sp_map.items()]))

    SeqIO.write(aln, outfile, "fasta-2line")

if __name__ == "__main__":
    main()
