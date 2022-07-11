from collections import Counter


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


def test_collect_exp_mut_freqs(coda):
    genome = "ATAGTCTAGCTGCATGACTGATCC"
    # genome = list("ATA GTC TAG CTG CAT GAC TGA TCC")
    #                  s   f     s f   s   s   s 
    nucl_freqs, cxt_freqs = coda.collect_exp_mut_freqs(list(genome))
    assert nucl_freqs["all"] == Counter(genome[1: -1])
    assert nucl_freqs["syn"] == {"A": 2, "C": 5, "G": 3, "T": 1}
    assert nucl_freqs["ff"] == {"C": 3, "G": 3}

    cxts = [genome[i: i+3] for i in range(len(genome)-2)]
    assert cxt_freqs["all"] == Counter(cxts)
    assert cxt_freqs["syn"] == {"TAG": 1, "TCT": 3, "GCT": 1, "TGC": 3, "ATG": 1, "ACT": 1, "GAT": 1}
    assert cxt_freqs["ff"] =={"TCT": 3, "TGC": 3}


def test_extract_mutations_simple():
    # TODO
    pass
