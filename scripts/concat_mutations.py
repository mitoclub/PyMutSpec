# USAGE: script.py [PATH_TO_REPLICS...] OUT_PATH_TO_MERGED_MUTATIONS

import re
import sys


def main():
    filepaths = sys.argv[1:-1]
    out = sys.argv[-1]
    # filepaths = ["tmp/evolve_sc/seqfile_sample-0000.fasta.mutations",
    #              "tmp/evolve_sc/seqfile_sample-0001.fasta.mutations"]
    # out = "tmp/evolve_sc/m.tsv"

    fout = open(out, "w")
    header_exists = False

    for i, fp in enumerate(filepaths):
        try:
            replica = int(re.search("sample-(\d{4}).fasta", fp).group(1))
        except:
            replica = i
        with open(fp) as fin:
            for line in fin:
                if line.startswith("Mut"):
                    header = line.strip().split("\t")
                    header.append("Replica")
                    if not header_exists:
                        fout.write("\t".join(header) + "\n")
                        header_exists = True
                else:
                    row = line.strip().split("\t")
                    row.append(str(replica))
                    fout.write("\t".join(row) + "\n")
    fout.close()


if __name__ == "__main__":
    main()
