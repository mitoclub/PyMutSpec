import os
import re
import glob
import sys

import click
import tqdm
import numpy as np
import pandas as pd


def read_marginal_probabilities(filepath, gene_name, model="F81"):
    nucleotides = list("ACGT")
    character = re.search(f"marginal_probabilities\.character_(.+)\.model_{model}\.tab", filepath)
    states = pd.read_csv(filepath, sep="\t")
    states.rename(columns={"node": "Node"}, inplace=True)
    states["State"] = [states.columns[i] for i in (np.argmax(states.iloc[:, 1:].values, axis=1) + 1)]
    for nucl in nucleotides:
        if nucl in states.columns:
            states["p_" + nucl] = states[nucl]
            states.drop(nucl, axis=1, inplace=True)
        else:
            states["p_" + nucl] = 0.0
    states["Site"] = int(character.groups()[0])
    states["Part"] = gene_name
    states = states[['Node', 'Part', 'Site', 'State', 'p_A', 'p_C', 'p_G', 'p_T']]
    return states


def extract_genename(indir):
    gn = indir.strip("/").split("/")[-1]
    gn = gn.removesuffix("_pastml")
    return gn


def read_gene_lens(dirpath):
    gene_lens = dict()
    for filepath in glob.glob(os.path.join(dirpath, "*")):
        gname = filepath.split("/")[-1].split(".")[0]
        with open(filepath) as fin:
            fin.readline()
            glen = len(fin.readline().strip())
            gene_lens[gname] = glen
    return gene_lens


def fill_gaps(states: pd.DataFrame, gene_lens):
    """
    Params
    ------
    states: pd.DataFrame
        all characters states for one gene
    """
    # return states

    assert states.Part.nunique() == 1
    gaps = []
    gene = states.iloc[0, 1]
    gl = gene_lens[gene]
    character_sites = set(states.Site.values)
    gap_sites = set(range(1, gl + 1)) - set(character_sites)
    for node in states.Node.unique():
        for site in gap_sites:
            gaps.append({
                "Node": node,
                "Part": gene,
                "Site": site,
                "State": "-",
                "p_A": 0.0,
                "p_C": 0.0,
                "p_G": 0.0,
                "p_T": 0.0,
            })
    gaps_df = pd.DataFrame(gaps)
    states_full = pd.concat([states, gaps_df])

    cur_gene_lens = states_full.groupby("Part").Site.max().to_dict()
    for gene in cur_gene_lens:
        if gene_lens[gene] != cur_gene_lens[gene]:
            raise RuntimeError(
                f"Len of gene are not equal to alignment length: "
                f"{gene_lens[gene]} != {cur_gene_lens[gene]}"
            )
    return states_full


def read_data(indir, gene_lens: dict):
    """
    Params
    ------
    indir: str
        path to dir, containing files with states for each separated character
    """
    gene = []
    gene_name = extract_genename(indir)
    for path in glob.glob(os.path.join(indir, "*")):
        if not "marginal_probabilities" in path:
            continue
        df = read_marginal_probabilities(path, gene_name)
        gene.append(df)
    gene_df = pd.concat(gene)
    gene_df_full = fill_gaps(gene_df, gene_lens)
    return gene_df_full


@click.command("formatter", help="reformat pastml output to usual states table")
@click.argument("dirs", nargs=-1, type=click.Path(True, False), required=True)
@click.option("--aln", "aln_dir", required=True, type=click.Path(True), help="path to directory with gene alignment files")
@click.option("--outpath", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(dirs, aln_dir, outpath):
    print(f"Input dirs:\n{dirs}")
    gene_lens = read_gene_lens(aln_dir)
    genome_states = []
    for data_dir in tqdm.tqdm(dirs, "Reading genes"):
        gene_states = read_data(data_dir, gene_lens)
        genome_states.append(gene_states)
    
    genome_df = pd.concat(genome_states)
    print("Sorting...", file=sys.stderr)
    genome_df: pd.DataFrame = genome_df.sort_values(["Node", "Part", "Site"])
    print("Writing...", file=sys.stderr)
    genome_df.to_csv(outpath, index=None, sep="\t", float_format="%.5f")


if __name__ == "__main__":
    main()
    # main(
    #     ["./data/pastml_n/ATP6_pastml",], 
    #     "./data/example_nematoda/alignments_nematoda_clean", 
    #     "./data/example_nematoda/genes_states.pastml.tsv"
    # )
