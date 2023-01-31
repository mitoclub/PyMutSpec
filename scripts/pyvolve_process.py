#!/usr/bin/env python3
from typing import List

import click
import pyvolve
import numpy as np
import pandas as pd
from Bio import SeqIO

from pymutspec.annotation import transcriptor
from pymutspec.io import read_rates

DEFAULT_REPLICS = 10
DEFAULT_GENCODE = 2


def get_rates(path_to_mutspec, eps=1e-3):
    ms = pd.read_csv(path_to_mutspec, sep="\t")
    ms["Mut"] = ms["Mut"].str.translate(transcriptor)
    ms["Mut"] = ms["Mut"].str.replace(">", "")
    ms["MutSpec"] = ms["MutSpec"] + eps
    rates = ms.set_index("Mut")["MutSpec"].to_dict()
    return rates


def get_root_seq(path_to_fasta):
    root_seq = None
    for rec in SeqIO.parse(path_to_fasta, format="fasta"):
        if rec.name.startswith("RN"):
            root_seq = str(rec.seq)
            break
    return root_seq


def codon_masking(seq, codons):
    """
    Arguments
    ---------
    seq: str
        root seq
    codons: array
        codon mask, 1-based
    """
    if len(seq) != codons:
        raise RuntimeError("lenghts of codon mask and seq are not equal")

    masked_codons = []
    codons = set(codons)
    j = 1
    for i in range(0, len(seq), 3):
        codon = seq[i: i+3]
        if len(codon) != 3:
            continue

        if j in codons:
            masked_codons.append(codon)
        j += 1

    return "".join(masked_codons)


def codon_unmasking(seq_base, codons, aln: List[str]):
    """
    Arguments
    ---------
    seq_base: str
        root seq
    codons: array
        codon mask, 1-based
    aln: List[str]
        simulated sequences without masked codons
    """
    if len(seq_base) != codons:
        raise RuntimeError("lenghts of codon mask and seq are not equal")

    unmasked_codons = []
    codons = set(codons)
    n_aln = len(aln)
    j = 1
    for i in range(0, len(seq_base), 3):
        codon_basic = seq_base[i: i+3]
        aln_codons = [x[i: i+3] for x in aln]
        if len(codon_basic) != 3:
            continue

        if j in codons:
            unmasked_codons.append(aln_codons)
        else:
            unmasked_codons.append([codon_basic for _ in range(n_aln)])
        j += 1

    for codon_pos in range(len(unmasked_codons[0])):
        
    return ["".join() for x in range(len(unmasked_codons[0]))]


@click.command("MutSel simulation", help="")
@click.option("-a", "--alignment", "path_to_mulal", type=click.Path(True), help="")
@click.option("-t", "--tree", "path_to_tree", type=click.Path(True), help="")
@click.option("-s", "--spectra", "path_to_mutspec", type=click.Path(True), help="")
@click.option("-o", "--out", type=click.Path(writable=True), help="Output sequences alignment (fasta)")
@click.option("--outcount", type=click.Path(writable=True), default=None, show_default=True, help="")
@click.option("-r", "--replics", "number_of_replics", default=DEFAULT_REPLICS, show_default=True, type=int, help="")
@click.option("-w", "--write_anc", is_flag=True, help="")
@click.option("-c", "--gencode", default=DEFAULT_GENCODE, show_default=True, help="")
@click.option("-l", "--scale_tree", default=1., show_default=True, help="")
@click.option("--rates", type=click.Path(True), default=None, help="path to rates from iqtree that will be used for positions masking")
def main(path_to_mulal, path_to_tree, path_to_mutspec, out, outcount, number_of_replics, write_anc, gencode, scale_tree, rates):
    tree = pyvolve.read_tree(file=path_to_tree, scale_tree=scale_tree)
    custom_mutation_asym = get_rates(path_to_mutspec)
    if rates is None:
        mask = codons = None
    else:
        mask = (read_rates(rates) > 1).astype(np.int8)  # gamma distribution category is one of [2,3,4,5]
        codons = np.where(np.reshape(mask[:len(mask) - len(mask)%3], (-1, 3)).sum(axis=1) > 0) + 1 # numbers of codons that changed
    codon_freqs = pyvolve.ReadFrequencies("codon", file=path_to_mulal, gencode=gencode, columns=codons).compute_frequencies(type="codon")
    model = pyvolve.Model("mutsel", {"state_freqs": codon_freqs, "mu": custom_mutation_asym}, gencode=gencode)
    root_seq = get_root_seq(path_to_mulal)
    partition = pyvolve.Partition(models=model, root_sequence=root_seq)
    evolver = pyvolve.Evolver(partitions=partition, tree=tree, gencode=gencode)

    for i in range(number_of_replics):
        print("Generating {} replica".format(i))
        evolver(
            seqfile=out.replace(".fasta", "_sample-{:04}.fasta".format(i)), 
            countfile=outcount,
            ratefile=None, infofile=None,
            write_anc=write_anc,
        )


if __name__ == "__main__":
    # main("-a ./tmp/evolve/alignment_checked.fasta -t ./tmp/evolve/iqtree_anc_tree.nwk -s ./tmp/ms12syn_.tsv -w -o ./tmp/evolve/seqfile.fasta".split())
    main()
