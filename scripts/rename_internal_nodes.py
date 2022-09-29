"""
script add names for internal nodes in newick file 
"""
import sys
from ete3 import PhyloTree


def main(path_to_tree, path_to_out):
    tree = PhyloTree(path_to_tree, format=1)

    i = 1
    for node in tree.traverse():
        if not node.name:
            node.name = "Node{}".format(i)
            i += 1

    tree.write(format=1, outfile=path_to_out)


if __name__ == "__main__":
    try:
        path_to_tree, path_to_out = sys.argv[1:]
        main(path_to_tree, path_to_out)
    except:
        print("ERROR\nUSAGE: script.py path_to_tree path_to_out_tree", file=sys.stderr)
