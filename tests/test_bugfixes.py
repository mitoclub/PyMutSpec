"""Regression tests for bugs fixed in 0.0.16."""
import os
import tempfile
import warnings

import numpy as np
import pandas as pd
import pytest
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pymutspec.annotation import (
    CodonAnnotation,
    calculate_mutspec,
    collapse_mutspec,
    complete_sbs192_columns,
    jackknife_spectra_sampling,
    calc_edgewise_spectra,
    get_tree_height,
    calc_phylocoefs,
)
from pymutspec.annotation.phylo_tree import Tree, TreeNode
from pymutspec.constants import possible_sbs192
from pymutspec.draw import plot_mutspec12
from pymutspec.io.states import GenomeStates
from pymutspec.io.gb import read_genbank_ref
from pymutspec.utils.logging import basic_logger


def test_prepare_codontable_invalid_raises():
    with pytest.raises(ValueError, match="not appropriate"):
        CodonAnnotation._prepare_codontable("not-a-table")


def test_get_syn_codons_returns_set(coda):
    result = coda.get_syn_codons("ATA", 1)
    assert isinstance(result, set)
    assert len(result) == 0


def test_cds_length_not_divisible_by_3(coda):
    seq3 = "ATGAAATAA"  # Met-Lys-Stop (9 nt); 3rd position of the stop is not synonymous
    seq4 = seq3 + "C"   # incomplete last codon
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        freqs3, _ = coda.collect_exp_mut_freqs(seq3, labels=["all", "syn"])
        freqs4, _ = coda.collect_exp_mut_freqs(seq4, labels=["all", "syn"])
    assert any(issubclass(w.category, UserWarning) for w in caught)
    # codon-aware (syn) counts ignore the extra base
    assert freqs3["syn"] == freqs4["syn"]
    # "all" can use the extra nucleotide as right-hand context of the last codon
    assert sum(freqs4["all"].values()) >= sum(freqs3["all"].values())


def test_collect_exp_muts_proba_missing_logger(coda):
    """CodonAnnotation has no logger; a bad mask must still raise ValueError."""
    genome = np.zeros((6, 4), dtype=float)
    genome[:, 0] = 1.0
    with pytest.raises(ValueError, match="same length"):
        coda.collect_exp_muts_proba(genome, phylocoef=1.0, mask=[1, 0])


def test_extract_mutations_simple(coda):
    g1 = np.array(list("ATGCTAGTA"))  # Met-Leu-Val
    g2 = np.array(list("ATGCTGGTA"))  # CTA -> CTG (Leu, fourfold)
    muts = coda.extract_mutations_simple(g1, g2)
    assert len(muts) == 1
    row = muts.iloc[0]
    assert row["Mut"] == "T[A>G]G"
    assert int(row["Label"]) == 2
    assert row["RefCodon"] == "CTA"
    assert row["AltCodon"] == "CTG"


def test_extract_mutations_simple_empty_has_columns(coda):
    g = np.array(list("ATGCTAGTA"))
    muts = coda.extract_mutations_simple(g, g)
    assert muts.empty
    assert "Mut" in muts.columns


def test_plot_mutspec_preserves_provided_axes():
    ms = pd.DataFrame({"Mut": ["C>A", "C>G", "C>T"], "MutSpec": [0.2, 0.3, 0.5]})
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plot_mutspec12(ms, ax=ax1, show=False, title="left")
    plot_mutspec12(ms, ax=ax2, show=False, title="right")
    assert fig.number in plt.get_fignums()
    assert ax1.get_title() == "left"
    assert ax2.get_title() == "right"
    plt.close(fig)


def test_complete_sbs192_columns_no_fragmentation():
    partial = pd.DataFrame({"A[C>A]A": [1.0], "T[G>T]T": [2.0]})
    complete = complete_sbs192_columns(partial)
    assert list(complete.columns) == possible_sbs192
    assert complete["A[C>A]A"].iloc[0] == 1.0
    assert complete["A[A>C]A"].iloc[0] == 0.0


def test_collapse_mutspec_accepts_obsnum_expnum():
    rows = []
    for sbs in possible_sbs192:
        rows.append({"Mut": sbs, "ObsNum": 1.0, "ExpNum": 1.0})
    ms192 = pd.DataFrame(rows)
    ms96 = collapse_mutspec(ms192)
    assert len(ms96) == 96
    assert "RawMutSpec" in ms96.columns


def test_jackknife_does_not_mutate_input():
    idx = pd.MultiIndex.from_tuples(
        [("R1", "A1"), ("R1", "A2"), ("R2", "A3")],
        names=["RefNode", "AltNode"],
    )
    obs = pd.DataFrame(1.0, index=idx, columns=possible_sbs192)
    exp = pd.DataFrame(1.0, index=pd.Index(["R1", "R2"], name="Node"), columns=possible_sbs192)
    obs_index_before = list(obs.index.names)
    exp_name_before = exp.index.name
    jackknife_spectra_sampling(obs, exp, frac=0.5, n=2)
    assert list(obs.index.names) == obs_index_before
    assert exp.index.name == exp_name_before


def test_nobs_cutoff_alias():
    idx = pd.MultiIndex.from_tuples(
        [("R1", "A1"), ("R1", "A2")],
        names=["RefNode", "AltNode"],
    )
    obs = pd.DataFrame(5.0, index=idx, columns=possible_sbs192)
    exp = pd.DataFrame(5.0, index=pd.Index(["R1"], name="Node"), columns=possible_sbs192)
    s1 = calc_edgewise_spectra(obs, exp, nmtypes_cutoff=0, nobs_cutoff=1, scale=False)
    s2 = calc_edgewise_spectra(obs, exp, nmtypes_cutoff=0, nobs_cuttof=1, scale=False)
    pd.testing.assert_frame_equal(s1, s2)


def test_tree_from_newick_string():
    t = Tree("(A:0.1,B:0.2)Root:0.0;")
    assert t.name in ("Root", "")
    leaves = sorted(n.name for n in t.iter_leaves())
    assert leaves == ["A", "B"]
    assert t.parent is None
    child = t.children[0]
    assert child.parent is t
    assert child.up is t
    edges = list(t.iter_edges())
    assert len(edges) == 2


def test_tree_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        Tree("definitely_not_a_tree_file.nwk")


def test_get_tree_height_zero_distances():
    root = TreeNode(name="R", dist=0.0)
    a = TreeNode(name="A", dist=0.0)
    b = TreeNode(name="B", dist=0.0)
    a._parent = b._parent = root
    root.children = [a, b]
    assert get_tree_height(root, "geom_mean") == 0.0
    assert get_tree_height(root, "mean") == 0.0
    coefs = calc_phylocoefs(root)
    assert coefs["R"] == 1.0


def test_genome_states_reads_states_not_gappy_file():
    states_txt = (
        "Node\tSite\tp_A\tp_C\tp_G\tp_T\n"
        "N1\t1\t1.0\t0.0\t0.0\t0.0\n"
        "N1\t2\t0.0\t1.0\t0.0\t0.0\n"
        "N1\t3\t0.0\t0.0\t1.0\t0.0\n"
    )
    gappy_txt = "99\n"
    with tempfile.TemporaryDirectory() as tmp:
        states_path = os.path.join(tmp, "states.tsv")
        gappy_path = os.path.join(tmp, "gappy.csv")
        with open(states_path, "w") as fh:
            fh.write(states_txt)
        with open(gappy_path, "w") as fh:
            fh.write(gappy_txt)
        gs = GenomeStates(states_path, path_to_gappy_sites=gappy_path, format="tsv")
        genome = gs.get_genome("N1")
        assert list(genome.index) == [1, 2, 3]
        assert gs.genome_size == 3


def test_genome_states_gappy_optional():
    states_txt = (
        "Node\tSite\tp_A\tp_C\tp_G\tp_T\n"
        "N1\t1\t1.0\t0.0\t0.0\t0.0\n"
        "N1\t2\t0.0\t1.0\t0.0\t0.0\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        states_path = os.path.join(tmp, "states.tsv")
        with open(states_path, "w") as fh:
            fh.write(states_txt)
        gs = GenomeStates(states_path, path_to_gappy_sites=None, format="tsv")
        assert gs.genome_size == 2


def test_basic_logger_does_not_duplicate_handlers():
    log1 = basic_logger()
    n1 = len(log1.handlers)
    log2 = basic_logger()
    assert log1 is log2
    assert len(log2.handlers) == n1


def test_read_genbank_ref_uses_gene_qualifier():
    gb = """\
LOCUS       TEST                 12 bp    DNA              UNK 01-JAN-1980
FEATURES             Location/Qualifiers
     source          1..12
                     /organism="test"
     CDS             1..12
                     /gene="cox1"
                     /codon_start=1
ORIGIN
        1 atgctagtaa tg
//
"""
    with tempfile.NamedTemporaryFile("w", suffix=".gb", delete=False) as fh:
        fh.write(gb)
        path = fh.name
    try:
        df = read_genbank_ref(path)
    finally:
        os.remove(path)
    assert "gene" in df.columns
    assert (df.loc[df["Type"] == "CDS", "gene"] == "cox1").all()


def test_filter_short_exp_seqs_drops_outlier():
    import sys
    sys.path.append("./scripts")
    from calculate_mutspec import filter_short_exp_seqs
    from pymutspec.constants import possible_sbs12, possible_sbs192

    cols = possible_sbs12 + possible_sbs192
    rows = []
    for node in [f"n{i}" for i in range(8)] + ["outlier"]:
        row = {c: 1.0 for c in cols}
        row["Node"] = node
        row["Label"] = "all"
        row["Gene"] = "g"
        if node == "outlier":
            for c in possible_sbs12:
                row[c] = 10000.0
        rows.append(row)
    flt = filter_short_exp_seqs(pd.DataFrame(rows))
    assert "outlier" not in set(flt.Node)
    assert len(flt) == 8
