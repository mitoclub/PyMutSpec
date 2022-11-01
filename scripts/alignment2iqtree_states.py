#!/usr/bin/env python3

"""
Leaves genomes must be written directly without partitioning - only positions in part!!
Only filenames need for this, but if there are preliminary aln concatenation it need additional step

Internal States must be rewritten to similar format 
"""

import os
import sys
from collections import defaultdict
from typing import Tuple

import click
from Bio import SeqIO

nucls = "ACGT"


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
        part = "1" if ngenes == 1 else os.path.basename(filepath).replace(".fna", "")
        fasta = SeqIO.parse(filepath, "fasta")
        for rec in fasta:
            node = rec.name
            history[node].append(part)
            seq = str(rec.seq)
            for site, state in enumerate(seq, 1):
                pos_data = [node, part, str(site), state]
                for nucl in nucls:
                    p = int(nucl == state)
                    pos_data.append(str(p))

                handle.write("\t".join(pos_data) + "\n")
        aln_lens[part] = len(seq)

    if ngenes > 1:
        # get max set of parts. Only for multiple files
        for node, parts in history.items():
            if len(parts) == ngenes:
                full_parts = parts.copy()
                break

        # fill missing genes by '-'. Only for multiple files
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
@click.argument("alignment", nargs=-1, type=click.Path(True))
@click.argument("states", nargs=1, type=click.Path(writable=True))
def main(alignment, states):
    parse_alignment_and_write_states(alignment, states)
    print(f"Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
    # main("./data/example_nematoda/alignments_nematoda_clean/", "/tmp/leaves_birds_states.tsv")
