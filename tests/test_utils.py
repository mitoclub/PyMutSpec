"""
Tests for utility functions in pymutspec.annotation and pymutspec.constants.

Covers:
- rev_comp / lbl2lbl_id / lbl_id2lbl (annotation.auxiliary)
- node_parent / iter_tree_edges (annotation.tree)
- get_iqr_bounds / filter_outlier_branches (annotation.spectra)
- complete_sbs192_columns / collapse_sbs192 (annotation.spectra)
- get_cossim / get_eucdist (annotation.spectra)
"""

import pytest
import pandas as pd

from pymutspec.annotation import (
    rev_comp, lbl2lbl_id, lbl_id2lbl,
    node_parent, iter_tree_edges,
    get_iqr_bounds, filter_outlier_branches,
    complete_sbs192_columns, collapse_sbs192,
    get_cossim, get_eucdist,
)
from pymutspec.constants import possible_sbs192, possible_sbs12


# ---------------------------------------------------------------------------
# rev_comp
# ---------------------------------------------------------------------------

class TestRevComp:
    def test_pyrimidine_mutation_unchanged_context(self):
        # C>A in context A_C = A[C>A]C
        # rev-comp: swap flanks (C,A) → C,A → complement → G,T
        # middle: [C>A] complement → [G>T]
        assert rev_comp("A[C>A]C") == "G[G>T]T"

    def test_involution(self):
        """Applying rev_comp twice returns the original string."""
        for sbs in ["A[C>A]C", "T[G>T]A", "C[A>G]T", "G[T>C]G"]:
            assert rev_comp(rev_comp(sbs)) == sbs

    def test_all_possible_sbs192_round_trip(self):
        """rev_comp is an involution on every element of possible_sbs192."""
        for sbs in possible_sbs192:
            assert rev_comp(rev_comp(sbs)) == sbs


# ---------------------------------------------------------------------------
# lbl_id2lbl / lbl2lbl_id
# ---------------------------------------------------------------------------

class TestLabelConversions:
    @pytest.mark.parametrize("lbl_id,lbl", [(0, "all"), (1, "syn"), (2, "ff")])
    def test_lbl_id2lbl(self, lbl_id, lbl):
        assert lbl_id2lbl(lbl_id) == lbl

    def test_lbl_id2lbl_invalid(self):
        with pytest.raises(NotImplementedError):
            lbl_id2lbl(99)

    @pytest.mark.parametrize("lbl,lbl_id", [
        ("all", 0), ("syn", 1), ("syn_c", 1), ("ff", 2), ("syn4f", 2),
    ])
    def test_lbl2lbl_id(self, lbl, lbl_id):
        assert lbl2lbl_id(lbl) == lbl_id

    def test_lbl2lbl_id_invalid(self):
        with pytest.raises(NotImplementedError):
            lbl2lbl_id("nonsyn")

    def test_round_trip(self):
        for lbl_id in (0, 1, 2):
            assert lbl2lbl_id(lbl_id2lbl(lbl_id)) == lbl_id


# ---------------------------------------------------------------------------
# node_parent / iter_tree_edges
# ---------------------------------------------------------------------------

class TestTreeUtils:
    def test_node_parent_root_returns_none(self, tree_rooted):
        assert node_parent(tree_rooted) is None

    def test_node_parent_leaf(self, tree_rooted):
        leaf = next(tree_rooted.iter_leaves())
        parent = node_parent(leaf)
        assert parent is not None
        assert leaf in parent.children

    def test_iter_tree_edges_count(self, tree_rooted):
        """Number of edges equals number of nodes minus one (tree property)."""
        n_nodes = sum(1 for _ in tree_rooted.traverse())
        edges = list(iter_tree_edges(tree_rooted))
        assert len(edges) == n_nodes - 1

    def test_iter_tree_edges_root_not_yielded_as_alt(self, tree_rooted):
        """The tree root must never appear as an alt_node."""
        root_name = tree_rooted.name
        for _, alt in iter_tree_edges(tree_rooted):
            assert alt.name != root_name

    def test_iter_tree_edges_ref_is_parent_of_alt(self, tree_rooted):
        for ref, alt in iter_tree_edges(tree_rooted):
            assert node_parent(alt) is ref


# ---------------------------------------------------------------------------
# get_iqr_bounds / filter_outlier_branches
# ---------------------------------------------------------------------------

class TestIQRAndOutlierFiltering:
    def test_get_iqr_bounds_simple(self):
        s = pd.Series([1, 2, 3, 4, 5, 100])
        lb, ub = get_iqr_bounds(s)
        assert lb < 1
        assert ub > 5
        assert 100 > ub  # 100 is an outlier

    def test_get_iqr_bounds_symmetric(self):
        s = pd.Series(range(1, 11))
        lb, ub = get_iqr_bounds(s)
        assert lb < 1
        assert ub > 10

    def test_filter_outlier_branches_removes_high_count(self):
        # Create a simple mutation table with one branch having many mutations
        normal = pd.DataFrame({
            "AltNode": ["A"] * 5 + ["B"] * 5 + ["C"] * 5,
            "Mut": ["A[C>A]T"] * 15,
            "ProbaMut": [1.0] * 15,
        })
        outlier = pd.DataFrame({
            "AltNode": ["OUTLIER"] * 100,
            "Mut": ["A[C>A]T"] * 100,
            "ProbaMut": [1.0] * 100,
        })
        obs_df = pd.concat([normal, outlier], ignore_index=True)
        filtered = filter_outlier_branches(obs_df, use_proba=True)
        assert "OUTLIER" not in filtered["AltNode"].values
        assert set(filtered["AltNode"].unique()).issubset({"A", "B", "C"})

    def test_filter_outlier_branches_no_proba(self):
        # Need enough normal branches for IQR to identify the outlier
        n_normal = 20
        normal = pd.DataFrame({
            "AltNode": [f"N{i}" for i in range(n_normal) for _ in range(5)],
            "Mut": ["A[C>A]T"] * (n_normal * 5),
            "ProbaMut": [1.0] * (n_normal * 5),
        })
        outlier = pd.DataFrame({
            "AltNode": ["OUTLIER"] * 1000,
            "Mut": ["A[C>A]T"] * 1000,
            "ProbaMut": [1.0] * 1000,
        })
        obs_df = pd.concat([normal, outlier], ignore_index=True)
        filtered = filter_outlier_branches(obs_df, use_proba=False)
        assert "OUTLIER" not in filtered["AltNode"].values


# ---------------------------------------------------------------------------
# complete_sbs192_columns / collapse_sbs192
# ---------------------------------------------------------------------------

class TestSbs192Utilities:
    def test_complete_sbs192_columns_fills_zeros(self):
        partial = pd.DataFrame(
            {"A[C>A]A": [1.0], "T[G>T]T": [2.0]},
        )
        complete = complete_sbs192_columns(partial)
        assert list(complete.columns) == possible_sbs192
        assert complete.shape == (1, 192)
        # Original values preserved
        assert complete["A[C>A]A"].iloc[0] == 1.0
        assert complete["T[G>T]T"].iloc[0] == 2.0
        # Missing columns filled with 0
        assert complete["A[A>C]A"].iloc[0] == 0.0

    def test_complete_sbs192_columns_already_complete(self):
        data = {sbs: [1.0] for sbs in possible_sbs192}
        df = pd.DataFrame(data)
        result = complete_sbs192_columns(df)
        assert list(result.columns) == possible_sbs192
        assert result.shape == (1, 192)

    def test_collapse_sbs192_to_12_shape(self):
        data = {sbs: [1.0] for sbs in possible_sbs192}
        df = pd.DataFrame(data)
        result = collapse_sbs192(df, to=12)
        assert list(result.columns) == possible_sbs12
        assert result.shape == (1, 12)

    def test_collapse_sbs192_to_12_sum_preserved(self):
        """Sum of 12-component spectrum equals sum of 192-component."""
        data = {sbs: [float(i)] for i, sbs in enumerate(possible_sbs192)}
        df = pd.DataFrame(data)
        result = collapse_sbs192(df, to=12)
        assert abs(result.values.sum() - df.values.sum()) < 1e-5

    def test_collapse_sbs192_invalid_to(self):
        data = {sbs: [1.0] for sbs in possible_sbs192}
        df = pd.DataFrame(data)
        with pytest.raises(NotImplementedError):
            collapse_sbs192(df, to=96)


# ---------------------------------------------------------------------------
# get_cossim / get_eucdist
# ---------------------------------------------------------------------------

class TestDistanceMetrics:
    def _make_dfs(self):
        cols = possible_sbs192
        a = pd.DataFrame(
            [[1.0] * 192, [0.5] * 192],
            index=["x", "y"],
            columns=cols,
        )
        b = pd.DataFrame(
            [[1.0] * 192, [1.0] * 192],
            index=["x", "z"],
            columns=cols,
        )
        return a, b

    def test_cossim_identical_vectors(self):
        cols = possible_sbs192
        df = pd.DataFrame([[1.0] * 192], index=["x"], columns=cols)
        result = get_cossim(df, df)
        assert abs(result["x"] - 1.0) < 1e-6

    def test_cossim_only_common_index(self):
        a, b = self._make_dfs()
        result = get_cossim(a, b)
        # Only "x" is common
        assert list(result.index) == ["x"]

    def test_cossim_empty_on_no_common_index(self):
        cols = possible_sbs192
        a = pd.DataFrame([[1.0] * 192], index=["x"], columns=cols)
        b = pd.DataFrame([[1.0] * 192], index=["y"], columns=cols)
        result = get_cossim(a, b)
        assert len(result) == 0

    def test_eucdist_same_vector_is_zero(self):
        cols = possible_sbs192
        df = pd.DataFrame([[0.5] * 192], index=["x"], columns=cols)
        result = get_eucdist(df, df)
        assert abs(result["x"]) < 1e-6

    def test_eucdist_only_common_index(self):
        a, b = self._make_dfs()
        result = get_eucdist(a, b)
        assert list(result.index) == ["x"]

    def test_eucdist_empty_on_no_common_index(self):
        cols = possible_sbs192
        a = pd.DataFrame([[1.0] * 192], index=["x"], columns=cols)
        b = pd.DataFrame([[1.0] * 192], index=["y"], columns=cols)
        result = get_eucdist(a, b)
        assert len(result) == 0
