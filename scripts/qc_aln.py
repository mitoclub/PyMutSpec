#!/usr/bin/env python3

import click
from Bio import SeqIO

@click.command("Fasta header coder", help="")
@click.option("-a", "--alignment", type=click.Path(True), help="Fasta nucleotide alignment")
@click.option("-o", "--outfile", type=click.Path(writable=True), help="Output fasta file with coded headers")
@click.option("-g", "--gap-frac",  default=0.3, type=float, help="Max fraction of gaps in sequence")
def main(alignment, outfile, gap_frac):
    aln = []
    for rec in SeqIO.parse(alignment, "fasta"):
        if rec.seq.count("-") / len(rec.seq) < gap_frac:
            aln.append(rec)

    SeqIO.write(aln, outfile, "fasta-2line")

if __name__ == "__main__":
    main()
