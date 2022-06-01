import os
import re
import glob
import sys

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
    gn = indir.split('/')[-1]
    gn = gn.removesuffix("_pastml")
    return gn


def read_data(indir):
    gene = []
    gene_name = extract_genename(indir)
    for path in glob.glob(os.path.join(indir, "*")):
        if not "marginal_probabilities" in path:
            continue
        df = read_marginal_probabilities(path, gene_name)
        gene.append(df)
    gene_df = pd.concat(gene)
    return gene_df


def read_gene_lens(dirpath):
    gene_lens = dict()
    for filepath in glob.glob(os.path.join(dirpath, "*")):
        gname = filepath.split("/")[-1].split(".")[0]
        with open(filepath) as fin:
            fin.readline()
            glen = len(fin.readline().strip())
            gene_lens[gname] = glen
    return gene_lens


def fill_gaps(states: pd.DataFrame, aln_dir, inplace=False):
    if not inplace:
        states = states.copy()
    gene_lens = read_gene_lens(aln_dir)
    gaps = []
    for gene in tqdm.tqdm(states.Part.unique(), "Filling genes gaps"):
        gl = gene_lens[gene]
        for node in states.Node.unique():
            ungapped_sites = states[(states.Node == node) & (states.Part == gene)].Site.values
            gapped_sites = set(range(1, gl + 1)) - set(ungapped_sites)
            for site in gapped_sites:
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
    states = pd.concat(states, gaps_df)
    cur_gene_lens = states.groupby("Part").Site.max().to_dict()
    assert gene_lens == cur_gene_lens
    return states


def main(dirs, outpath, aln_dir):
    genome_states = []
    for data_dir in tqdm.tqdm(dirs, "Reading genes"):
        gene_states = read_data(data_dir)
        genome_states.append(gene_states)
    
    genome_df = pd.concat(genome_states)
    fill_gaps(genome_df, aln_dir, inplace=True)
    print("Sorting...", file=sys.stderr)
    genome_df = genome_df.sort_values(["Node", "Part", "Site"])
    genome_df.to_csv(outpath, index=None, sep="\t")
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main(["./data/pastml_n/ATP6_pastml",], "/tmp/states.tsv", "./data/example_nematoda/alignments_nematoda_clean")
