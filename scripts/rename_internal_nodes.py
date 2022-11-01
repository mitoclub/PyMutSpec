#!/usr/bin/env python3

import sys
from ete3 import PhyloTree, PhyloNode

dist_formatter = "%0.8f"


def node_parent(node: PhyloNode):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def main(path_to_dist_tree, path_to_named_tree, path_to_out):
    """Iteratively going from leaves to root and add names to dist tree"""
    tree_dist = PhyloTree(path_to_dist_tree, format=0)
    tree_named = PhyloTree(path_to_named_tree, format=8)

    nd = len(tree_dist.get_cached_content())
    nn = len(tree_named.get_cached_content())
    tree_named.traverse()
    if nd != nn:
        raise RuntimeError(f"trees are different: {nd} and {nn} nodes")

    for leaf_dist in tree_dist.iter_leaves():
        node_dist = leaf_dist
        node_named = next(tree_named.iter_search_nodes(name=node_dist.name))

        while node_dist != tree_dist:
            pa_dist = node_parent(node_dist)
            pa_named = node_parent(node_named)
            if pa_dist.name and pa_dist.name != pa_named.name:
                raise RuntimeError("pa_dist.name != pa_named.name")
            pa_dist.name = pa_named.name

            node_dist = pa_dist
            node_named = pa_named

    nwk = tree_dist.write(format=1, outfile=None, dist_formatter=dist_formatter)
    nwk = nwk.replace(";", "ROOT;")  # add ROOT label, that cannot be added with PhyloTree
    with open(path_to_out, "w") as fout:
        fout.write(nwk)


if __name__ == "__main__":
    try:
        path_to_dist_tree, path_to_named_tree, path_to_out = sys.argv[1:]
        main(path_to_dist_tree, path_to_named_tree, path_to_out)
    except:
        print("ERROR\nUSAGE: script.py path_to_dist_tree path_to_named_tree path_to_out_tree", file=sys.stderr)
