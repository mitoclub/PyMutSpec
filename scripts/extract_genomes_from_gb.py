import random
import sys

from Bio import SeqIO
import tqdm

N = 300
p = 0.1
print(N, p)


def main():
    try:
        path_in = sys.argv[1]
        path_out = sys.argv[2]
    except:
        print("USAGE: script.py GB_FILE OUT_FASTA")

    records = SeqIO.parse(path_in, format="gb")
    nd1_seqs = []
    cytb_seqs = []
    for rec in tqdm.tqdm(records):
        if random.random() < p:
            for fea in rec.features:
                if fea.type == "CDS":
                    if fea.qualifiers["gene"][0] == "CYTB":
                        cytb_seqs.append(fea.extract(rec))
                    if fea.qualifiers["gene"][0] == "ND1":
                        nd1_seqs.append(fea.extract(rec))

        if len(nd1_seqs) == N:
            break
    
    with open(path_out.replace(".fasta", "_cytb.fasta"), "w") as fout:
        SeqIO.write(cytb_seqs, fout, "fasta")
    
    with open(path_out.replace(".fasta", "_nd1.fasta"), "w") as fout:
        SeqIO.write(nd1_seqs, fout, "fasta")


if __name__ == "__main__":
    main()
