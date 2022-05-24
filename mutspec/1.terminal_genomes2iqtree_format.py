"""
Leaves genomes must be written directly without partitioning - only positions in gene!!
Only filenames need for this, but if there are preliminary aln concatenation it need additional step

Internal States must be rewritten to similar format 
"""

import os
import sys
from collections import defaultdict
from typing import Tuple

import click
import tqdm
from Bio import SeqIO

from mutspec.utils import load_scheme, get_aln_files

NGENES = 12


def parse_alignment(files: list, scheme: dict, aln_dir, out) -> Tuple[str, int]:
    """
    read fasta files from scheme with alignments and write states to table

    return states table and full alignment length
    """
    handle = open(out, "w")
    columns = "Node Part Site State p_A p_C p_G p_T".split()
    handle.write("\t".join(columns) + "\n")
    aln_lens = []
    files = set(files)
    history = defaultdict(list)
    for part, gene_fn in tqdm.tqdm(scheme.items(), "Parts"):
        filepath = os.path.join(aln_dir, gene_fn)
        assert filepath in files, f"cannot find file {filepath} from scheme"
        fasta = SeqIO.parse(filepath, "fasta")
        for rec in fasta:
            node = rec.name
            history[node].append(part)
            seq = str(rec.seq)
            for site, state in enumerate(seq, 1):
                pos_data = [node, str(part), str(site), state]
                for nucl in "ACGT":
                    p = int(nucl == state)
                    pos_data.append(str(p))

                handle.write("\t".join(pos_data) + "\n")
        aln_lens.append(len(seq))

    # get max set of parts
    for node, parts in history.items():
        if len(parts) == NGENES:
            full_parts = parts.copy()
            break
    
    # fill missing genes by '-'
    for node, parts in history.items():
        if len(parts) != NGENES:
            unseen_parts = set(full_parts).difference(parts)
            for unp in unseen_parts:
                print(f"Gap filling for node {node}, part {unp}...", file=sys.stderr)
                for site in range(1, aln_lens[unp - 1] + 1):
                    pos_data = [node, str(unp), str(site), "-", "0", "0", "0", "0"]
                    handle.write("\t".join(pos_data) + "\n")

    handle.close()
    full_aln_len = sum(aln_lens)
    return full_aln_len


@click.command("formatter", help="reformat alignment to states table")
@click.option("--aln", "aln_dir", required=True, type=click.Path(True), help="path to directory with gene alignment files")
@click.option("--scheme", "scheme_path", required=True, type=click.Path(True), help="path to scheme that contain gene splitting info of alignment")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(aln_dir, scheme_path, out):
    aln_files = get_aln_files(aln_dir)
    scheme = load_scheme(scheme_path)
    assert len(scheme) > 0
    print(scheme)
    aln_len = parse_alignment(aln_files, scheme, aln_dir, out)


if __name__ == "__main__":
    main()
    # main("data/interim/alignments_birds_clean_clean", "data/interim/scheme_birds_genes.nex", "data/interim/leaves_birds_states.tsv")
