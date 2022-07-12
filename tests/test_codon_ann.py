from collections import Counter
import pytest

def test_get_syn_codons(coda):
    assert len(coda.get_syn_codons("ATA", 1)) == 0
    assert len(coda.get_syn_codons("ATA", 2)) == 1
    assert len(coda.get_syn_codons("CTA", 0)) == 1
    assert len(coda.get_syn_codons("CTA", 2)) == 3
    assert len(coda.get_syn_codons("TAA", 2)) == 0
    assert len(coda.get_syn_codons("AAA", 2)) == 1
    assert len(coda.get_syn_codons("ACC", 2)) == 3
    assert len(coda.get_syn_codons("CGG", 0)) == 0
    assert len(coda.get_syn_codons("CGG", 2)) == 3


def test_get_mut_type(coda):
    with pytest.raises(ValueError):
        coda.get_mut_type("TAT", "TAC", 1)
    with pytest.raises(ValueError):
        coda.get_mut_type("TAT", "TAT", 1)
    with pytest.raises(ValueError):
        coda.get_mut_type("TAT", "TAC", 3)
    
    assert coda.get_mut_type("TAT", "TAC", 2) == (1, "Y", "Y")
    assert coda.get_mut_type("CAA", "CAC", 2) == (0, "Q", "H")
    assert coda.get_mut_type("CTT", "CTA", 2) == (2, "L", "L")
    assert coda.get_mut_type("CTG", "TTG", 0) == (1, "L", "L")
    assert coda.get_mut_type("CCT", "ACT", 0) == (0, "P", "T")
    assert coda.get_mut_type("TAC", "TAA", 2) == (-1, "Y", "*")
    assert coda.get_mut_type("AGA", "AGC", 2) == (-2, "*", "S")
    assert coda.get_mut_type("AGA", "AGG", 2) == (-3, "*", "*")


def test_collect_exp_mut_freqs(coda):
    genome = "ATAGTCTAGCTGCATGACTGATCC"
    # genome = list("ATA GTC TAG CTG CAT GAC TGA TCC")
    #                  s   f     s f   s   s   s 
    exp_sbs12_freqs, exp_sbs192_freqs = coda.collect_exp_mut_freqs(genome)
    expected_sbs12_all = dict()
    for nuc1, freq in Counter(genome[1: -1]).items():
        for nuc2 in "ACGT":
            if nuc1 != nuc2:
                expected_sbs12_all[f"{nuc1}>{nuc2}"] = freq
    assert exp_sbs12_freqs["all"] == expected_sbs12_all
    assert exp_sbs12_freqs["ff"] == {"C>A": 1, "C>G": 1, "C>T": 1, "G>A": 1, "G>C": 1, "G>T": 1}
    assert exp_sbs12_freqs["syn"] == {"A>G": 2, "C>A": 1, "C>G": 1, "C>T": 3, "G>A": 1, "G>C": 1, "G>T": 1, "T>C": 1}

    cxts = [genome[i: i+3] for i in range(len(genome)-2)]
    expected_sbs192_all = dict()
    for cxt, freq in Counter(cxts).items():
        nuc1 = cxt[1]
        for nuc2 in "ACGT":
            if nuc1 != nuc2:
                expected_sbs192_all[f"{cxt[0]}[{nuc1}>{nuc2}]{cxt[-1]}"] = freq
    assert exp_sbs192_freqs["all"] == expected_sbs192_all
    assert exp_sbs192_freqs["syn"] == {
        "T[A>G]G": 1, "G[C>T]T": 1,
        "T[C>A]T": 1, "T[C>G]T": 1, "T[C>T]T": 1, 
        "T[G>A]C": 1, "T[G>C]C": 1, "T[G>T]C": 1, 
        "A[C>T]T": 1, "A[T>C]G": 1, "G[A>G]T": 1
    }
    assert exp_sbs192_freqs["ff"] == {
        "T[C>A]T": 1, "T[C>G]T": 1, "T[C>T]T": 1, 
        "T[G>A]C": 1, "T[G>C]C": 1, "T[G>T]C": 1
    }


def test_extract_mutations_simple():
    # TODO
    pass
