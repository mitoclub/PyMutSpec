#!/usr/bin/env python3

"""
Leaves genomes must be written directly without partitioning - only positions in gene!!
Only filenames need for this, but if there are preliminary aln concatenation it need additional step

Internal States must be rewritten to similar format 
"""

import os
import sys
import glob
from collections import defaultdict
from typing import Tuple

import click
from Bio import SeqIO


def get_aln_files(path: str):
    assert os.path.isdir(path), "path is not directory"
    files = set(glob.glob(os.path.join(path, "*.fna")))
    return files


def parse_alignment_and_write_states(files: list, outfile) -> Tuple[str, int]:
    """
    read fasta files from scheme with alignments and write states to table

    return states table and full alignment length
    """
    print(f"Processing...", file=sys.stderr)
    handle = open(outfile, "w")
    ngenes = len(files)
    columns = "Node Part Site State p_A p_C p_G p_T".split()
    handle.write("\t".join(columns) + "\n")
    aln_lens = dict()
    files = set(files)
    history = defaultdict(list)
    for filepath in files:
        gene = os.path.basename(filepath).replace(".fna", "")
        fasta = SeqIO.parse(filepath, "fasta")
        for rec in fasta:
            node = rec.name
            history[node].append(gene)
            seq = str(rec.seq)
            for site, state in enumerate(seq, 1):
                pos_data = [node, str(gene), str(site), state]
                for nucl in "ACGT":
                    p = int(nucl == state)
                    pos_data.append(str(p))

                handle.write("\t".join(pos_data) + "\n")
        aln_lens[gene] = len(seq)

    # get max set of parts
    for node, parts in history.items():
        if len(parts) == ngenes:
            full_parts = parts.copy()
            break

    # fill missing genes by '-'
    for node, parts in history.items():
        if len(parts) != ngenes:
            unseen_parts = set(full_parts).difference(parts)
            for unp in unseen_parts:
                print(f"Gap filling for node {node}, part {unp}...", file=sys.stderr)
                for site in range(1, aln_lens[unp] + 1):
                    pos_data = [node, unp, str(site), "-", "0", "0", "0", "0"]
                    handle.write("\t".join(pos_data) + "\n")
    handle.close()


@click.command("formatter", help="reformat alignment to states table")
@click.option("--aln", required=True, type=click.Path(True), help="path to directory with gene alignment files")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(aln, out):
    aln_files = get_aln_files(aln)
    parse_alignment_and_write_states(aln_files, out)


if __name__ == "__main__":
    main()
    # main("./data/example_nematoda/alignments_nematoda_clean/", "/tmp/leaves_birds_states.tsv")
