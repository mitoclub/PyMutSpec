from queue import Queue
from statistics import geometric_mean

import numpy as np

# TODO merge with phylo_tree.py and remove this file; fix tests after these changes

def node_parent(node):
    """
    TODO remove this function and use node._parent directly AND rename node._parent to node.parent

    Return the parent node of *node*, or ``None`` if *node* is the root.

    Arguments
    ---------
    node
        A node in a phylogenetic tree.

    Return
    ------
    parent or None
        The immediate ancestor of *node*, or ``None`` when *node* has no
        ancestors (i.e. it is the root).
    """
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def iter_tree_edges(tree):
    """
    TODO integrate to Tree class

    Iterate over all directed edges (parent → child) in the tree via BFS.

    The root node itself is skipped; every other node produces exactly one
    ``(ref_node, alt_node)`` pair where *ref_node* is the parent and
    *alt_node* is the child.

    Arguments
    ---------
    tree
        Rooted phylogenetic tree.

    Yields
    ------
    ref_node
        Parent (reference) node of the edge.
    alt_node
        Child (alternative) node of the edge.
    """
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


def get_tree_height(tree, mode='geom_mean'):
    """
    Return the characteristic length of a (sub)tree as the distance from
    the root to its leaves.

    Arguments
    ---------
    tree
        Rooted phylogenetic tree or subtree.  Must not be named ``'ROOT'``.
    mode: str
        Aggregation method over leaf distances.  One of:

        - ``'mean'``      – arithmetic mean of leaf distances
        - ``'geom_mean'`` – geometric mean of leaf distances (default)
        - ``'max'``       – distance to the farthest leaf

    Return
    ------
    tree_len: float
        Characteristic length of the tree.

    Raises
    ------
    TypeError
        If *mode* is not one of the accepted values.
    """
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
        raise TypeError("mode must be 'mean', 'geom_mean' or 'max'")

    return md


def get_ingroup_root(tree):
    """
    Return the ingroup root of a binary rooted tree that contains an outgroup.

    The function assumes the tree root has exactly two children, one of which
    is a leaf (the outgroup).  If no leaf child is found the tree root itself
    is returned.

    Arguments
    ---------
    tree
        Rooted binary tree with an outgroup leaf attached to the root.

    Return
    ------
    ingrp
        Root of the ingroup clade.

    Raises
    ------
    AssertionError
        If the tree root does not have exactly two children.
    """
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


def calc_phylocoefs(tree):
    """
    Calculate a phylogenetic coefficient for every node in the tree.

    The coefficient for a node is ``1 - d / tree_len``, where *d* is the
    distance from the node to its closest leaf and *tree_len* is the
    geometric-mean leaf distance of the ingroup root.  Values are capped so
    that the minimum coefficient is > 0 (i.e. ``d / tree_len`` is capped
    at 0.99999).

    Arguments
    ---------
    tree
        Rooted binary phylogenetic tree (with an outgroup leaf, optional).

    Return
    ------
    phylocoefs: dict[str, float]
        Mapping from node name to its phylogenetic coefficient.
    """
    ingroup = get_ingroup_root(tree)
    tree_height = get_tree_height(ingroup, 'geom_mean')
    root_phylocoef = 1 - min(0.999, ingroup.get_closest_leaf()[1] / tree_height)
    phylocoefs = {ingroup.name: root_phylocoef}
    for node in ingroup.iter_descendants():
        _closest, d = node.get_closest_leaf()
        phylocoefs[node.name] = 1 - min(0.99999, d / tree_height)
    return phylocoefs
