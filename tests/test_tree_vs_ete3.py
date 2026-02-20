"""
Tests comparing the custom Tree/TreeNode implementation against ete3's
PhyloTree/PhyloNode using the same newick file.
"""
import pytest
from ete3 import PhyloTree as Ete3PhyloTree

from pymutspec.annotation.phylo_tree import Tree
from pymutspec.annotation import iter_tree_edges

PATH = "./tests/data/treefile_rooted.nwk"


@pytest.fixture(scope="module")
def ete3_tree():
    return Ete3PhyloTree(PATH, format=1)


@pytest.fixture(scope="module")
def custom_tree():
    return Tree(PATH, format=1)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def test_root_name(ete3_tree, custom_tree):
    """Root node names must be identical."""
    assert custom_tree.name == ete3_tree.name


def test_root_dist(ete3_tree, custom_tree):
    """Root branch lengths must be identical."""
    assert round(custom_tree.dist, 8) == round(ete3_tree.dist, 8)


def test_leaf_count(ete3_tree, custom_tree):
    """len(tree) returns the number of leaves in both implementations."""
    assert len(custom_tree) == len(ete3_tree)


def test_total_node_count(ete3_tree, custom_tree):
    """get_cached_content() must return the same number of nodes."""
    assert len(custom_tree.get_cached_content()) == len(ete3_tree.get_cached_content())


# ---------------------------------------------------------------------------
# Node naming
# ---------------------------------------------------------------------------

def test_all_leaf_names_match(ete3_tree, custom_tree):
    """The sorted set of leaf names must be identical."""
    ete3_leaves = sorted(n.name for n in ete3_tree.iter_leaves())
    custom_leaves = sorted(n.name for n in custom_tree.iter_leaves())
    assert custom_leaves == ete3_leaves


def test_all_node_names_match(ete3_tree, custom_tree):
    """The sorted set of all node names must be identical."""
    ete3_nodes = sorted(n.name for n in ete3_tree.traverse())
    custom_nodes = sorted(n.name for n in custom_tree.traverse())
    assert custom_nodes == ete3_nodes


def test_search_nodes_by_name(ete3_tree, custom_tree):
    """search_nodes(name=...) must return a node with the correct name."""
    node_name = "Node1"
    ete3_match = ete3_tree.search_nodes(name=node_name)
    custom_match = custom_tree.search_nodes(name=node_name)
    assert len(custom_match) == len(ete3_match) == 1
    assert custom_match[0].name == ete3_match[0].name


# ---------------------------------------------------------------------------
# Branch / edge iteration
# ---------------------------------------------------------------------------

def test_iter_tree_edges_count(ete3_tree, custom_tree):
    """iter_tree_edges must yield the same number of edges for both trees."""
    ete3_edges = list(iter_tree_edges(ete3_tree))
    custom_edges = list(iter_tree_edges(custom_tree))
    assert len(custom_edges) == len(ete3_edges)


def test_iter_tree_edges_node_names(ete3_tree, custom_tree):
    """The sorted (ref, alt) name pairs from iter_tree_edges must be identical."""
    def edge_name_pairs(tree):
        return sorted((ref.name, alt.name) for ref, alt in iter_tree_edges(tree))

    assert edge_name_pairs(custom_tree) == edge_name_pairs(ete3_tree)


def test_iter_descendants_count(ete3_tree, custom_tree):
    """iter_descendants must yield the same number of nodes."""
    ete3_desc = list(ete3_tree.iter_descendants())
    custom_desc = list(custom_tree.iter_descendants())
    assert len(custom_desc) == len(ete3_desc)


def test_iter_leaves_names(ete3_tree, custom_tree):
    """iter_leaves must yield leaves with the same names (sorted)."""
    ete3_leaves = sorted(n.name for n in ete3_tree.iter_leaves())
    custom_leaves = sorted(n.name for n in custom_tree.iter_leaves())
    assert custom_leaves == ete3_leaves


# ---------------------------------------------------------------------------
# Branch lengths
# ---------------------------------------------------------------------------

def test_all_branch_lengths_match(ete3_tree, custom_tree):
    """
    For every named node, the branch length reported by the custom tree must
    match ete3's value (to 8 decimal places).
    """
    ete3_dists = {n.name: n.dist for n in ete3_tree.traverse() if n.name}
    custom_dists = {n.name: n.dist for n in custom_tree.traverse() if n.name}

    assert set(custom_dists.keys()) == set(ete3_dists.keys())
    for name in ete3_dists:
        assert round(custom_dists[name], 8) == round(ete3_dists[name], 8), (
            f"Branch length mismatch for node '{name}': "
            f"custom={custom_dists[name]}, ete3={ete3_dists[name]}"
        )
