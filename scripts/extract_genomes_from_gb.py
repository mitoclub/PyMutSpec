from collections import defaultdict
import random
import sys

from Bio import SeqIO
import tqdm

N = 300
p = 1
print(N, p)

# genes = ["CYTB", "ND1"]
genes = []

def main():
    try:
        path_in = sys.argv[1]
        path_out = sys.argv[2]
    except:
        print("USAGE: script.py GB_FILE OUT_FASTA")

    records = SeqIO.parse(path_in, format="gb")
    seqs = defaultdict(list)

    for rec in tqdm.tqdm(records):
        if random.random() < p:
            for fea in rec.features:
                if fea.type == "CDS":
                    for g in genes:
                        if fea.qualifiers["gene"][0] == g:
                            seqs[g].append(fea.extract(rec))
                elif fea.type == "D-loop":
                    new_rec = fea.extract(rec)
                    new_rec.id = rec.id
                    new_rec.name = rec.description
                    new_rec.description = rec.description
                    seqs["D-loop"].append(new_rec)

        if len(genes) and len(seqs[genes[0]]) == N:
            break
    
    for g in genes + ["D-loop"]:
        with open(path_out.replace(".fasta", f"_{g}.fasta"), "w") as fout:
            SeqIO.write(seqs[g], fout, "fasta")



if __name__ == "__main__":
    main()
