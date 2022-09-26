#!/usr/bin/env python3

# Reformat ancestral states of RAxML to IQTREE style
# USAGE: script.py PATH_TO_RAxML_STATES PATH_TO_OUT_IQTREE_STATES

import sys

nucls = list("ACGT")
header = ["Node", "Site", "State", "p_A", "p_C", "p_G", "p_T"]


def agrmax(a):
    xmax = 0
    imax = None
    for i, x in enumerate(a):
        if x > 0.5:
            xmax = x
            imax = i
            break
        if x > xmax:
            xmax = x
            imax = i
    return imax


def main(path_to_raxml_states, out_path):
    with open(out_path, "w") as fout:
        fout.write("\t".join(header) + "\n")
        with open(path_to_raxml_states) as fin:
            node = fin.readline().strip()
            i = 1
            for line in fin:
                line = line.strip()
                if not line:
                    node = fin.readline().strip()
                    i = 1
                else:
                    probas = line.split()
                    most_probable_state = nucls[agrmax(map(float, probas))]
                    row = [node, str(i), most_probable_state] + probas
                    fout.write("\t".join(row) + "\n")
                    i += 1


if __name__ == "__main__":
    # path_to_raxml_states = "data/RAxML_marginalAncestralProbabilities.tsv"
    # out_path = "data/out.state"
    try:
        main(sys.argv[1], sys.argv[2])
    except:
        RuntimeError("USAGE: script.py PATH_TO_RAxML_STATES PATH_TO_OUT_IQTREE_STATES")
