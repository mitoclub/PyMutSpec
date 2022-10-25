import pyvolve
import pandas as pd

path_to_mulal_fasta = "../Downloads/alignment_checked.fasta"
path_to_tree = "./tmp/dist.nwk"
path_to_mutspec = "./tmp/ms12syn_.tsv"

ms = pd.read_csv(path_to_mutspec, sep="\t")
ms["Mut"] = ms["Mut"].str.replace(">", "")
ms["MutSpec"] = ms["MutSpec"] + 0.01
custom_mutation_asym = ms.set_index("Mut")["MutSpec"].to_dict()

my_tree = pyvolve.read_tree(file=path_to_tree)
f = pyvolve.ReadFrequencies("codon", file=path_to_mulal_fasta)
codon_freqs = f.compute_frequencies(type="codon")
params = {"state_freqs": codon_freqs, "mu": custom_mutation_asym}
my_model = pyvolve.Model("mutsel", params)
my_partition = pyvolve.Partition(models=my_model, size=200)
my_evolver = pyvolve.Evolver(partitions=my_partition, tree=my_tree)
my_evolver(ratefile = "custom_ratefile.txt", infofile = "custom_infofile.txt", seqfile
= "custom_seqfile.fasta" )

