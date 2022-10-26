from email.policy import default
import click
import pyvolve
import pandas as pd
from Bio import SeqIO

DEFAULT_REPLICS = 10
DEFAULT_ALGO = "exponentiation"


def get_rates(path_to_mutspec):
    ms = pd.read_csv(path_to_mutspec, sep="\t")
    ms["Mut"] = ms["Mut"].str.replace(">", "")
    ms["MutSpec"] = ms["MutSpec"] + 0.01
    rates = ms.set_index("Mut")["MutSpec"].to_dict()
    return rates


def get_root_seq(path_to_fasta):
    root_seq = None
    for rec in SeqIO.parse(path_to_fasta, format="fasta"):
        if rec.name.startswith("RN"):
            root_seq = str(rec.seq)#[:144]
            break
    return root_seq


@click.command("MutSel simulation", help="")
@click.option("-a", "--alignment", "path_to_mulal", type=click.Path(True), help="")
@click.option("-t", "--tree", "path_to_tree", type=click.Path(True), help="")
@click.option("-s", "--spectra", "path_to_mutspec", type=click.Path(True), help="")
@click.option("-r", "--replics", "number_of_replics", default=DEFAULT_REPLICS, show_default=True, type=int, help="")
@click.option("-l", "--algorithm", default=DEFAULT_ALGO, show_default=True, type=click.Choice(["exponentiation", "Gillespie"]), help="")
@click.option("-w", "--write_anc", is_flag=True, help="")
def main(path_to_mulal, path_to_tree, path_to_mutspec, number_of_replics, algorithm, write_anc):
    algorithm = 0 if algorithm == "exponentiation" else 1
    tree = pyvolve.read_tree(file=path_to_tree)
    custom_mutation_asym = get_rates(path_to_mutspec)
    codon_freqs = pyvolve.ReadFrequencies("codon", file=path_to_mulal).compute_frequencies(type="codon")
    model = pyvolve.Model("mutsel", {"state_freqs": codon_freqs, "mu": custom_mutation_asym})
    root_seq = get_root_seq(path_to_mulal)
    partition = pyvolve.Partition(models=model, root_sequence=root_seq)
    evolver = pyvolve.Evolver(partitions=partition, tree=tree)

    # for _ in range(number_of_replics):
    evolver(
        seqfile="tmp/evolve/custom_seqfile.fasta",
        countfile="tmp/evolve/countfile.txt",
        ratefile=None, infofile=None,
        write_anc=write_anc,
        # algorithm=algorithm,
        algorithm=0,
    )


if __name__ == "__main__":
    main("-a ./tmp/alignment_checked.fasta -t ./tmp/dist.nwk -s ./tmp/ms12syn_.tsv -w".split())
