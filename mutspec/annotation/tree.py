from queue import Queue

import numpy as np
from ete3 import PhyloTree


def node_parent(node):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def iter_tree_edges(tree: PhyloTree):
    discovered_nodes = set()
    discovered_nodes.add(tree.name)
    Q = Queue()
    Q.put(tree)

    while not Q.empty():
        cur_node = Q.get()
        for child in cur_node.children:
            Q.put(child)

        if cur_node.name not in discovered_nodes:
            discovered_nodes.add(cur_node.name)
            alt_node = cur_node
            ref_node = node_parent(alt_node)
            yield ref_node, alt_node


def get_farthest_leaf(tree: PhyloTree, quantile=None):
    """
    TODO check if tree is rooted
    """
    if quantile is None:
        _, md = tree.get_farthest_leaf()
        return md
    elif isinstance(quantile, (float, int)) and 0 <= quantile <= 1:
        distances = []
        for leaf in tree.iter_leaves():
            d = tree.get_distance(leaf)
            distances.append(d)
        md = np.quantile(distances, quantile)
    else:
        raise TypeError(f"quantile must be int, float or None, got {type(quantile)}")
    return md
