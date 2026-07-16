from pymutspec.annotation import get_ingroup_root, get_tree_height, calc_phylocoefs


def test_get_ingroup_root(tree_rooted):
    ingroup = get_ingroup_root(tree_rooted)
    assert ingroup.name == 'Node1'
    assert round(ingroup.dist, 4) == 0.0739


def test_get_tree_height(tree_rooted):
    ingroup = get_ingroup_root(tree_rooted)

    h1 = get_tree_height(tree_rooted, 'geom_mean')
    h2 = get_tree_height(ingroup, 'geom_mean')

    assert round(h1, 4) == 0.3435
    assert round(h2, 4) == 0.2689


def test_calc_phylocoefs(tree_rooted):
    phylocoefs = calc_phylocoefs(tree_rooted)

    ingroup = get_ingroup_root(tree_rooted)
    tl = get_tree_height(ingroup, 'geom_mean')
    
    node = 'Node1'
    n1_d = tree_rooted.search_nodes(name=node)[0].get_closest_leaf()[1]

    assert round(phylocoefs[node], 5) == round(1 - n1_d / tl, 5)
