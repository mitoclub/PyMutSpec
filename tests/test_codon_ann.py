from collections import Counter


def test_get_syn_number(coda):
    assert coda.get_syn_number("ATA", 1) == 0
    assert coda.get_syn_number("ATA", 2) == 1
    assert coda.get_syn_number("CTA", 0) == 1
    assert coda.get_syn_number("CTA", 2) == 3
    assert coda.get_syn_number("TAA", 2) == 0
    assert coda.get_syn_number("AAA", 2) == 1
    assert coda.get_syn_number("ACC", 2) == 3
    assert coda.get_syn_number("CGG", 0) == 0
    assert coda.get_syn_number("CGG", 2) == 3


def test_collect_obs_mut_freqs(coda):
    genome = "ATAGTCTAGCTGCATGACTGATCC"
    # genome = list("ATA GTC TAG CTG CAT GAC TGA TCC")
    #                  s   f     s f   s   s   s 
    nucl_freqs, cxt_freqs = coda.collect_obs_mut_freqs(list(genome))
    assert nucl_freqs["all"] == Counter(genome[1: -1])
    assert nucl_freqs["syn"] == {"A": 2, "C": 5, "G": 3, "T": 1}
    assert nucl_freqs["ff"] == {"C": 1, "G": 1}

    cxts = [genome[i: i+3] for i in range(len(genome)-2)]
    assert cxt_freqs["all"] == Counter(cxts)
    assert cxt_freqs["syn"] == {"TAG": 1, "TCT": 3, "GCT": 1, "TGC": 3, "ATG": 1, "ACT": 1, "GAT": 1}
    assert cxt_freqs["ff"] =={"TCT": 1, "TGC": 1}


def test_extract_mutations_simple():
    # TODO
    pass
