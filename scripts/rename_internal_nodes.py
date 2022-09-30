import sys
from ete3 import PhyloTree, PhyloNode


def node_parent(node: PhyloNode):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def main(path_to_dist_tree, path_to_named_tree, path_to_out):
    # Iteratively going from leaves to root and add names to dist tree
    tree_dist = PhyloTree(path_to_dist_tree, format=0)
    tree_named = PhyloTree(path_to_named_tree, format=8)
    
    nd = len(tree_dist.get_cached_content())
    nn = len(tree_named.get_cached_content())
    tree_named.traverse()
    if nd != nn:
        raise RuntimeError(f"trees are different: {nd} and {nn} nodes")

    for leaf_dist in tree_dist.iter_leaves():
        leaf_named = tree_named.search_nodes(name=leaf_dist.name)
        pa_named = node_parent(leaf_named)
        pa_dist = node_parent(leaf_dist)
        pa_dist.name = pa_named.name

        




if __name__ == "__main__":
    try:
        # path_to_dist_tree, path_to_named_tree, path_to_out = sys.argv[1:]
        path_to_dist_tree, path_to_named_tree, path_to_out = "tmp/dist.nwk", "tmp/named.nwk", "tmp/out.nwk"
        main(path_to_dist_tree, path_to_named_tree, path_to_out)
    except:
        print("ERROR\nUSAGE: script.py path_to_dist_tree path_to_named_tree path_to_out_tree", file=sys.stderr)
