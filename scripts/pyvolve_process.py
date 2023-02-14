#!/usr/bin/env python3
from typing import Dict

import click
import pyvolve
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    if len(seq) // 3 < max(codons):
        raise ValueError("codon mask don't fit to sequence")

    masked_codons = []
    codons = set(codons)
    codon_idx = 1
    for i in range(0, len(seq), 3):
        if codon_idx in codons:
            codon = seq[i: i+3]
            assert len(codon) == 3
            masked_codons.append(codon)
        codon_idx += 1

    return "".join(masked_codons)


def codon_unmasking(seq_common, codons, aln_dct: Dict[str, str]) -> Dict[str, str]:
    """
    Arguments
    ---------
    seq_common: str
        root seq
    codons: array
        codon mask, 1-based
    aln: List[str]
        simulated sequences without masked codons
    """
    if len(seq_common) // 3 < max(codons):
        raise ValueError("codon mask don't fit to sequence")

    unmasked_codons = []
    codons = set(codons)
    aln = list(aln_dct.values())
    n_aln = len(aln)
    codon_idx_in_common = 1
    shift = 0
    for i in range(0, len(seq_common), 3):
        codon_basic = seq_common[i: i+3]
        if len(codon_basic) != 3:
            raise ValueError("Codon length != 3")

        if codon_idx_in_common in codons:
            aln_codons = [x[i - shift: i - shift + 3] for x in aln]
            unmasked_codons.append(aln_codons)
        else:
            unmasked_codons.append([codon_basic for _ in range(n_aln)])
            shift += 3
        codon_idx_in_common += 1
    
    seqs = ["".join([x[i] for x in unmasked_codons]) for i in range(len(unmasked_codons[0]))]
    headers = list(aln_dct.keys())
    unmasked_dct = dict(zip(headers, seqs))
    return unmasked_dct


def write_seqs(seqdict: Dict[str, str], seqfile, seqfmt="fasta-2line"):
    alignment = []
    for entry in seqdict:
        seq_entry = Seq( seqdict[entry])
        seq_object = SeqRecord(seq_entry, id = entry, description = "", annotations={"molecule_type": "DNA"}) 
        alignment.append(seq_object)     

    SeqIO.write(alignment, seqfile, seqfmt)


@click.command("MutSel simulation", help="")
@click.option("-a", "--alignment", "path_to_mulal", type=click.Path(True), help="")
@click.option("-t", "--tree", "path_to_tree", type=click.Path(True), help="")
@click.option("-s", "--spectra", "path_to_mutspec", type=click.Path(True), help="")
@click.option("--rates", type=click.Path(True), default=None, help="path to rates from iqtree that will be used for positions masking")
@click.option("-o", "--out", type=click.Path(writable=True), help="Output sequences alignment (fasta)")
@click.option("--outcount", type=click.Path(writable=True), default=None, help="")
@click.option("-r", "--replics", "nreplics", default=DEFAULT_REPLICS, show_default=True, type=int, help="")
@click.option("-w", "--write_anc", is_flag=True, help="")
@click.option("-c", "--gencode", default=DEFAULT_GENCODE, show_default=True, help="")
@click.option("-l", "--scale_tree", default=1., show_default=True, help="")
def main(path_to_mulal, path_to_tree, path_to_mutspec, out, outcount, nreplics, write_anc, gencode, scale_tree, rates):
    if rates is None:
        mask = codons = columns = None
        scale_tree_factor = scale_tree
    else:
        mask = (read_rates(rates) > 1).astype(np.int8)  # gamma distribution category is one of [2,3,4,5]
        codons = np.where(np.reshape(mask[:len(mask) - len(mask)%3], (-1, 3)).sum(axis=1) > 0)[0] + 1 # indexes of codons that changed
        columns=list(codons)
        scale_tree_factor = scale_tree * len(mask) / (len(codons) * 3)
        print("Tree scaling factor = {:.2f}".format(scale_tree_factor))

    tree = pyvolve.read_tree(file=path_to_tree, scale_tree=scale_tree_factor)
    custom_mutation_asym = get_rates(path_to_mutspec)
    codon_freqs = pyvolve.ReadFrequencies("codon", file=path_to_mulal, gencode=gencode, columns=columns).compute_frequencies(type="codon")
    model = pyvolve.Model("mutsel", {"state_freqs": codon_freqs, "mu": custom_mutation_asym}, gencode=gencode)
    root_seq = get_root_seq(path_to_mulal)
    root_seq_masked = root_seq if codons is None else codon_masking(root_seq, codons)
    partition = pyvolve.Partition(models=model, root_sequence=root_seq_masked)
    evolver = pyvolve.Evolver(partitions=partition, tree=tree, gencode=gencode)

    for i in range(nreplics):
        print("Generating {} replica".format(i))
        seqfile = out.replace(".fasta", "_sample-{:04}.fasta".format(i))
        if codons is None:
            evolver(
                seqfile=seqfile, 
                countfile=outcount,
                ratefile=None, infofile=None,
                write_anc=write_anc,
            )
        else:
            evolver(
                seqfile=None,
                countfile=outcount,
                ratefile=None, infofile=None,
                write_anc=write_anc,
            )
            aln_dct = evolver.get_sequences(True)
            unmasked_aln = codon_unmasking(root_seq, codons, aln_dct)
            write_seqs(unmasked_aln, seqfile)


if __name__ == "__main__":
    # main("-a data/exposure/human_cytb/pyvolve/mulal.fasta.clean -t data/exposure/human_cytb/pyvolve/tree.nwk.ingroup " 
    #      "-s data/exposure/human_cytb/ms/ms12syn.tsv -o data/tmp/seqfile.fasta -w -r 2 -c 2 "
    #      "--rates data/exposure/human_cytb/CYTB.rate".split())
    main()
