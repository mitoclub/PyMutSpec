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


def get_ingroup_root(tree: PhyloTree, outgrp='OUTGRP'):
    found_outgrp, found_ingroup = False, False
    ingroup = None
    for n in tree.children:
        if n.name == outgrp:
            found_outgrp = True
        else:
            found_ingroup = True
            ingroup = n
    if found_ingroup and found_outgrp:
        return ingroup
    else:
        raise Exception('Cannot extract ingroup root')


def calc_phylocoefs(tree: PhyloTree, outgrp='OUTGRP'):
    tree_len = get_tree_len(get_ingroup_root(tree, outgrp ), 'geom_mean')
    phylocoefs = {}
    for node in tree.iter_descendants():
        d = node.get_closest_leaf()[1]
        phylocoefs[node.name] = 1 - min(0.9999, d / tree_len)
    return phylocoefs


def get_tree_outgrp_name(tree: PhyloTree):
    c1, c2 = tree.children
    done = False
    if len(c1.children) == 0:
        done = True
        outgrp = c1
    if not done and len(c2.children) == 0:
        done = True
        outgrp = c2
    if not done:
        raise Exception("Cannot get outgroup from tree")
        
    return outgrp.name
