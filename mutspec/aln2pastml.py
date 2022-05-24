"""
Leaves genomes must be written directly without partitioning - only positions in gene!!
Only filenames need for this, but if there are preliminary aln concatenation it need additional step

Internal States must be rewritten to similar format 

TODO rewrite input - only aln-files(!!!) BUT WE'll LOOSE THE PARTS INDEXES

"""

import os
from typing import Tuple

import click
import tqdm
import pandas as pd
from Bio import SeqIO

from mutspec.utils import load_scheme, get_aln_files

NGENES = 12


def parse_alignment(files: list, scheme: dict, aln_dir) -> Tuple[str, int]:
    """
    read fasta files from scheme with alignments and write states to table

    return states table and full alignment length
    """
    files = set(files)
    data = []
    for part, gene_fn in tqdm.tqdm(scheme.items(), "Parts"):
        filepath = os.path.join(aln_dir, gene_fn)
        assert filepath in files, f"cannot find file {filepath} from scheme"
        fasta = SeqIO.parse(filepath, "fasta")
        for rec in fasta:
            node = rec.name
            seq = str(rec.seq)
            for site, state in enumerate(seq, 1):
                pos_data = [node, part, site, state]
                data.append(pos_data)

    df = pd.DataFrame(data, columns="Node Part Site State".split())
    print("Sorting...")
    df = df.sort_values(["Node", "Part", "Site"])
    
    print("Positioning...")
    genes_length = df.groupby("Part").Site.max().to_dict()
    part_pos = {}
    for part in genes_length:
        pos = 0
        for i in range(1, part):
            pos += genes_length[i]
        part_pos[part] = pos

    df["PosInAln"] = df.Part.map(part_pos)
    df["PosInAln"] = df["PosInAln"] + df["Site"]
    df.drop(["Part", "Site"], axis=1, inplace=True)
    
    print("Pivoting...")
    pidf = df.pivot("Node", "PosInAln", "State")
    return pidf


@click.command("formatter", help="reformat alignment to states table")
@click.option("--aln", "aln_dir", required=True, type=click.Path(True), help="path to directory with gene alignment files")
@click.option("--scheme", "scheme_path", required=True, type=click.Path(True), help="path to scheme that contain gene splitting info of alignment")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(aln_dir, scheme_path, out):
    aln_files = get_aln_files(aln_dir)
    scheme = load_scheme(scheme_path)
    print(scheme)
    df = parse_alignment(aln_files, scheme, aln_dir)
    df.to_csv(out, sep="\t")


if __name__ == "__main__":
    main()
    # main("data/example_birds/aln", "data/example_birds/scheme_birds_genes.nex", "data/leaves_birds.pastml.tsv")
