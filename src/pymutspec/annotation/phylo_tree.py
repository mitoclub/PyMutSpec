"""
Custom phylogenetic tree classes replacing ete3 dependency.
Uses BioPython for newick format parsing.
"""
from io import StringIO

from Bio import Phylo as _BioPhylo


class TreeNode:
    """
    A phylogenetic tree node with an ete3-compatible interface.
    Instances form a tree through parent/child references.
    """

    def __init__(self, name="", dist=0.0):
        self.name = name if name is not None else ""
        self.dist = dist if dist is not None else 0.0
        self.children = []
        self._parent = None

    # ------------------------------------------------------------------
    # Tree navigation
    # ------------------------------------------------------------------

    def is_leaf(self):
        return len(self.children) == 0

    def traverse(self):
        """Yield all nodes (self first, then descendants)."""
        yield self
        for child in self.children:
            yield from child.traverse()

    def iter_descendants(self):
        """Yield all descendant nodes (not self)."""
        for child in self.children:
            yield child
            yield from child.iter_descendants()

    def iter_ancestors(self):
        """Yield ancestor nodes from parent to root."""
        node = self._parent
        while node is not None:
            yield node
            node = node._parent

    def iter_leaves(self):
        """Yield all leaf nodes reachable from self."""
        if self.is_leaf():
            yield self
        else:
            for child in self.children:
                yield from child.iter_leaves()

    # ------------------------------------------------------------------
    # Distance computation
    # ------------------------------------------------------------------

    def get_distance(self, target):
        """Return cumulative branch length from self to *target* (descendant)."""
        def _find(node, target, acc):
            if node is target:
                return acc
            for child in node.children:
                result = _find(child, target, acc + child.dist)
                if result is not None:
                    return result
            return None

        d = _find(self, target, 0.0)
        if d is None:
            raise ValueError(f"Node '{target.name}' is not a descendant of '{self.name}'")
        return d

    def get_farthest_leaf(self):
        """Return *(leaf, distance)* for the leaf farthest from self."""
        best_leaf, best_dist = None, -1.0
        for leaf in self.iter_leaves():
            d = self.get_distance(leaf)
            if d > best_dist:
                best_dist = d
                best_leaf = leaf
        return best_leaf, best_dist

    def get_closest_leaf(self):
        """Return *(leaf, distance)* for the leaf closest to self."""
        best_leaf, best_dist = None, float("inf")
        for leaf in self.iter_leaves():
            d = self.get_distance(leaf)
            if d < best_dist:
                best_dist = d
                best_leaf = leaf
        return best_leaf, best_dist

    # ------------------------------------------------------------------
    # Node lookup
    # ------------------------------------------------------------------

    def get_cached_content(self):
        """Return a dict keyed by every node reachable from self (ete3 compat)."""
        return {node: None for node in self.traverse()}

    def search_nodes(self, name=None, **kwargs):
        """Return list of all nodes whose attributes match the given criteria."""
        results = []
        for node in self.traverse():
            if name is not None and node.name != name:
                continue
            if all(getattr(node, k, None) == v for k, v in kwargs.items()):
                results.append(node)
        return results

    def iter_search_nodes(self, name=None, **kwargs):
        """Yield all nodes whose attributes match the given criteria."""
        for node in self.traverse():
            if name is not None and node.name != name:
                continue
            if all(getattr(node, k, None) == v for k, v in kwargs.items()):
                yield node

    # ------------------------------------------------------------------
    # Newick serialisation
    # ------------------------------------------------------------------

    def _to_newick(self, dist_formatter=None):
        fmt = dist_formatter if dist_formatter else "%g"
        if self.is_leaf():
            return f"{self.name}:{fmt % self.dist}"
        children_str = ",".join(c._to_newick(dist_formatter) for c in self.children)
        return f"({children_str}){self.name}:{fmt % self.dist}"

    def write(self, format=1, outfile=None, dist_formatter=None):  # noqa: A002
        """
        Serialise tree to newick.

        Parameters
        ----------
        format : int
            Newick format hint (kept for API compatibility with ete3; currently
            ignored – the tree is always written in a standard newick format
            that includes all node names and branch lengths).
        outfile : str or None
            If given, also write the string to this file path.
        dist_formatter : str or None
            printf-style format string for branch lengths (e.g. ``"%0.8f"``).

        Returns
        -------
        str
            Newick string.
        """
        children_str = ",".join(c._to_newick(dist_formatter) for c in self.children)
        fmt = dist_formatter if dist_formatter else "%g"
        nwk = f"({children_str}){self.name}:{fmt % self.dist};"
        if outfile:
            with open(outfile, "w") as fh:
                fh.write(nwk)
        return nwk

    def __len__(self):
        """Return the number of leaf nodes (mirrors ete3 behaviour)."""
        return sum(1 for _ in self.iter_leaves())

    def __repr__(self):
        return f"TreeNode(name={self.name!r}, dist={self.dist}, children={len(self.children)})"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _bio_clade_to_node(clade):
    """Recursively convert a BioPython Clade into a TreeNode tree."""
    node = TreeNode(
        name=clade.name if clade.name is not None else "",
        dist=clade.branch_length if clade.branch_length is not None else 0.0,
    )
    for child_clade in clade.clades:
        child_node = _bio_clade_to_node(child_clade)
        child_node._parent = node
        node.children.append(child_node)
    return node


# ---------------------------------------------------------------------------
# Public constructor – mirrors ete3's PhyloTree(path, format=N)
# ---------------------------------------------------------------------------

class Tree(TreeNode):
    """
    Load a phylogenetic tree from a newick file.

    Parameters
    ----------
    newick_path : str
        Path to a newick-format tree file.
    format : int
        Newick format hint (kept for API compatibility with ete3; this
        parameter is currently ignored – BioPython's newick parser handles
        all common newick variants automatically).
    """

    def __init__(self, newick_path, format=1):  # noqa: A002
        with open(newick_path) as fh:
            tree_str = fh.read().strip()

        bio_tree = _BioPhylo.read(StringIO(tree_str), "newick")
        root = _bio_clade_to_node(bio_tree.root)

        # Initialise self as the root node (TreeNode.__init__ not called via
        # super because we copy attributes from the parsed root directly).
        self.name = root.name
        self.dist = root.dist
        self.children = root.children
        self._parent = None
        # Re-point children's _parent to self (they currently point to root).
        for child in self.children:
            child._parent = self

