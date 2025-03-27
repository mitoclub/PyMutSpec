from queue import Queue
from statistics import geometric_mean

import numpy as np
from ete3 import PhyloTree, PhyloNode


def node_parent(node: PhyloNode):
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


def get_tree_len(tree: PhyloTree, mode='geom_mean'):
    '''
    TODO check if tree is rooted 

    Params:
        - mode: str - calculate 'mean', 'geom_mean' or 'max' of distribution of len from current node to leaves
    '''
    assert tree.name != 'ROOT'

    if mode == 'max':
        _, md = tree.get_farthest_leaf()
    elif mode in ['mean', 'geom_mean']:
        distances_to_leaves = []
        for leaf in tree.iter_leaves():
            d = tree.get_distance(leaf)
            distances_to_leaves.append(d)
        
        if mode == 'mean':
            md = np.mean(distances_to_leaves)
        elif mode == 'geom_mean':
            md = geometric_mean(distances_to_leaves)

    else:
        raise TypeError(f"mode must be 'mean', 'geom_mean' or 'max'")

    return md


def get_ingroup_root(tree: PhyloTree) -> PhyloTree:
    assert len(tree.children) == 2, 'Tree must be binary'
    found_outgroup = False
    for node in tree.children:
        if node.is_leaf():
            found_outgroup = True
        else:
            ingrp = node

    if found_outgroup:
        return ingrp
    else:
        return tree


def calc_phylocoefs(tree: PhyloTree):
    tree_len = get_tree_len(get_ingroup_root(tree), 'geom_mean')
    phylocoefs = {tree.name: 1 - min(0.999, tree.get_closest_leaf()[1] / tree_len)}
    for node in tree.iter_descendants():
        _closest, d = node.get_closest_leaf()
        phylocoefs[node.name] = 1 - min(0.99999, d / tree_len)
    return phylocoefs
