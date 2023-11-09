import pytest
from pymutspec.annotation import get_tree_outgrp_name, get_ingroup_root, get_tree_len, calc_phylocoefs


def test_get_tree_outgrp_name(tree_rooted):
    outgrp = get_tree_outgrp_name(tree_rooted)
    assert outgrp == 'Acanthisitta_chloris'


def test_get_ingroup_root(tree_rooted):
    ingroup = get_ingroup_root(tree_rooted, 'Acanthisitta_chloris')
    assert ingroup.name == 'Node1'
    assert round(ingroup.dist, 4) == 0.0739


def test_get_tree_len(tree_rooted):
    ingroup = get_ingroup_root(tree_rooted, 'Acanthisitta_chloris')

    l1 = get_tree_len(tree_rooted, 'geom_mean')
    l2 = get_tree_len(ingroup, 'geom_mean')

    assert round(l1, 4) == 0.3435
    assert round(l2, 4) == 0.2689


def test_calc_phylocoefs(tree_rooted):
    phylocoefs = calc_phylocoefs(tree_rooted, 'Acanthisitta_chloris')

    ingroup = get_ingroup_root(tree_rooted, 'Acanthisitta_chloris')
    tl = get_tree_len(ingroup, 'geom_mean')
    
    node = 'Node1'
    n1_d = tree_rooted.search_nodes(name=node)[0].get_closest_leaf()[1]

    assert round(phylocoefs[node], 5) == round(1 - n1_d / tl, 5)
