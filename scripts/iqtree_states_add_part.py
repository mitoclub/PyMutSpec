#!/usr/bin/env python3

"""
script add "Part" column to states (tsv) file
"""
import sys


def main(path_to_states, path_to_out):
    with open(path_to_out, "w") as fout:
        with open(path_to_states) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                row = line.strip().split()
                if line.startswith("Node\t"):
                    header = [row[0]] + ["Part"] + row[1:]
                    fout.write("\t".join(header) + "\n")
                else:
                    row = [row[0]] + ["1"] + row[1:]
                    fout.write("\t".join(row) + "\n")


if __name__ == "__main__":
    try:
        path_to_states, path_to_out = sys.argv[1:]
        main(path_to_states, path_to_out)
    except:
        print("ERROR\nUSAGE: script.py path_to_states path_to_out_states", file=sys.stderr)
