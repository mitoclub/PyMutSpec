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
    labels = ['all', 'syn', 'ff', 'nonsyn']
    exp_sbs12_freqs, exp_sbs192_freqs = coda.collect_exp_mut_freqs(genome, labels=labels)
    expected_sbs12_all = dict()
    for nuc1, freq in Counter(genome[1: -1]).items():
        for nuc2 in "ACGT":
            if nuc1 != nuc2:
                expected_sbs12_all[f"{nuc1}>{nuc2}"] = freq
    expected_sbs12_nonsyn = {
        "A>C": 5, "A>G": 3, "A>T": 5,
        "C>A": 4, "C>G": 4, "C>T": 2,
        "G>A": 4, "G>C": 4, "G>T": 4,
        "T>A": 7, "T>C": 6, "T>G": 7,
    }
    expected_sbs12_syn = {"A>G": 2, "C>A": 1, "C>G": 1, "C>T": 3, "G>A": 1, "G>C": 1, "G>T": 1, "T>C": 1}
    expected_sbs12_ff = {"C>A": 1, "C>G": 1, "C>T": 1, "G>A": 1, "G>C": 1, "G>T": 1}
    diff_of_all_and_syn = {sbs12: x-expected_sbs12_syn.get(sbs12, 0) for sbs12,x in expected_sbs12_all.items()}
    assert exp_sbs12_freqs["all"] == expected_sbs12_all
    assert exp_sbs12_freqs["syn"] == expected_sbs12_syn
    assert exp_sbs12_freqs["ff"] == expected_sbs12_ff
    assert exp_sbs12_freqs["nonsyn"] == expected_sbs12_nonsyn
    assert exp_sbs12_freqs["nonsyn"] == diff_of_all_and_syn


    cxts = [genome[i: i+3] for i in range(len(genome)-2)]
    expected_sbs192_all = dict()
    for cxt, freq in Counter(cxts).items():
        nuc1 = cxt[1]
        for nuc2 in "ACGT":
            if nuc1 != nuc2:
                expected_sbs192_all[f"{cxt[0]}[{nuc1}>{nuc2}]{cxt[-1]}"] = freq
    expected_sbs192_syn = {
        "T[A>G]G": 1, "G[C>T]T": 1,
        "T[C>A]T": 1, "T[C>G]T": 1, "T[C>T]T": 1, 
        "T[G>A]C": 1, "T[G>C]C": 1, "T[G>T]C": 1, 
        "A[C>T]T": 1, "A[T>C]G": 1, "G[A>G]T": 1
    }
    expected_sbs192_ff = {
        "T[C>A]T": 1, "T[C>G]T": 1, "T[C>T]T": 1, 
        "T[G>A]C": 1, "T[G>C]C": 1, "T[G>T]C": 1
    }
    expected_sbs192_nonsyn = {s: x-expected_sbs192_syn.get(s, 0) for s,x in expected_sbs192_all.items() if x-expected_sbs192_syn.get(s, 0) > 0}
    assert exp_sbs192_freqs["all"] == expected_sbs192_all
    assert exp_sbs192_freqs["syn"] == expected_sbs192_syn
    assert exp_sbs192_freqs["ff"]  == expected_sbs192_ff
    assert exp_sbs192_freqs["nonsyn"]  == expected_sbs192_nonsyn


def test_collect_exp_mut_freqs_on_real_gene(coda, states_most_probable):
    genome = states_most_probable.get_genome('Node1')['1']
    labels = ['all', 'syn', 'nonsyn']
    exp_sbs12_freqs, exp_sbs192_freqs = coda.collect_exp_mut_freqs(genome, labels=labels)
    expected_sbs12_all = dict()
    for nuc1, freq in Counter(genome[1: -1]).items():
        for nuc2 in "ACGT":
            if nuc1 != nuc2:
                expected_sbs12_all[f"{nuc1}>{nuc2}"] = freq

    cxts = ["".join(genome[i: i+3]) for i in range(len(genome)-2)]
    expected_sbs192_all = dict()
    for cxt, freq in Counter(cxts).items():
        nuc1 = cxt[1]
        for nuc2 in "ACGT":
            if nuc1 != nuc2:
                expected_sbs192_all[f"{cxt[0]}[{nuc1}>{nuc2}]{cxt[-1]}"] = freq

    expected_sbs12_nonsyn = {sbs12: x - exp_sbs12_freqs['syn'].get(sbs12, 0) \
                             for sbs12,x in expected_sbs12_all.items()}
    expected_sbs192_nonsyn = {sbs192: x - exp_sbs192_freqs['syn'].get(sbs192, 0) \
                              for sbs192,x in expected_sbs192_all.items() \
                                if x - exp_sbs192_freqs['syn'].get(sbs192, 0) > 0}
    
    assert exp_sbs12_freqs["all"] == expected_sbs12_all
    assert exp_sbs12_freqs["nonsyn"] == expected_sbs12_nonsyn
    assert exp_sbs192_freqs["all"] == expected_sbs192_all
    assert exp_sbs192_freqs["nonsyn"]  == expected_sbs192_nonsyn


def test_collect_exp_mut_freqs_on_real_gene_proba(coda, states):
    genome = states.get_genome('Node1')['1']
    labels = ['all', 'syn', 'nonsyn']
    exp_sbs12_freqs, exp_sbs192_freqs = coda.collect_exp_mut_freqs_proba(
        genome, phylocoef=1., labels=labels)

    expected_sbs12_nonsyn = {sbs12: x - exp_sbs12_freqs['syn'].get(sbs12, 0) \
                             for sbs12, x in exp_sbs12_freqs['all'].items()}
    expected_sbs192_nonsyn = {sbs192: x - exp_sbs192_freqs['syn'].get(sbs192, 0) \
                              for sbs192, x in exp_sbs192_freqs['all'].items() \
                                if x - exp_sbs192_freqs['syn'].get(sbs192, 0) > 0}
    
    assert exp_sbs12_freqs["nonsyn"] == expected_sbs12_nonsyn
    assert exp_sbs192_freqs["nonsyn"]  == expected_sbs192_nonsyn


def test_collect_exp_muts_on_real_gene(coda, states_most_probable):
    genome = states_most_probable.get_genome('Node1')['1']
    labels = ['all', 'syn', 'ff', 'nonsyn']
    exp192 = coda.collect_exp_muts(genome, labels=labels)
    _, exp_sbs192_freqs = coda.collect_exp_mut_freqs(genome, labels=labels)
    
    exp_freqs_all = exp192[exp192.Label == 'all'].groupby('Mut').Mut.count().to_dict()
    exp_freqs_syn = exp192[exp192.Label == 'syn'].groupby('Mut').Mut.count().to_dict()
    exp_freqs_ff = exp192[exp192.Label == 'syn4f'].groupby('Mut').Mut.count().to_dict()
    exp_freqs_nonsyn = exp192[exp192.Label == 'nonsyn'].groupby('Mut').Mut.count().to_dict()

    assert exp_sbs192_freqs['all'] == exp_freqs_all
    assert exp_sbs192_freqs['syn'] == exp_freqs_syn
    assert exp_sbs192_freqs['ff'] == exp_freqs_ff
    assert exp_sbs192_freqs['nonsyn'] == exp_freqs_nonsyn


def test_collect_exp_muts_on_real_gene_proba(coda, states):
    genome = states.get_genome('Node1')['1']
    labels = ['all', 'syn', 'ff', 'nonsyn']
    exp192 = coda.collect_exp_muts_proba(
        genome, phylocoef=1., labels=labels)
    _, exp_sbs192_freqs = coda.collect_exp_mut_freqs_proba(
        genome, phylocoef=1., labels=labels)
    
    exp_freqs_all = exp192[exp192.Label == 'all'].groupby('Mut').Proba.sum().to_dict()
    exp_freqs_syn = exp192[exp192.Label == 'syn'].groupby('Mut').Proba.sum().to_dict()
    exp_freqs_ff = exp192[exp192.Label == 'syn4f'].groupby('Mut').Proba.sum().to_dict()
    exp_freqs_nonsyn = exp192[exp192.Label == 'nonsyn'].groupby('Mut').Proba.sum().to_dict()

    assert exp_sbs192_freqs['all'] == exp_freqs_all
    assert exp_sbs192_freqs['syn'] == exp_freqs_syn
    assert exp_sbs192_freqs['ff'] == exp_freqs_ff
    assert exp_sbs192_freqs['nonsyn'] == exp_freqs_nonsyn


def test_extract_mutations_simple():
    # TODO
    pass

# use to test collect_exp_muts_proba and collect_exp_mut_freqs_proba (results must be comparable and equal)
# np.all(coda.collect_exp_muts_proba(states.get_genome("Node4705")["1"], 1, mut_proba_cutoff=0.05).groupby(["Label", "Mut"]).Proba.sum().unstack()[possible_sbs192].fillna(0) == \
#     pd.DataFrame(coda.collect_exp_mut_freqs_proba(states.get_genome("Node4705")["1"], 1, mut_proba_cutoff=0.05)[1]).T[possible_sbs192].rename(index={"ff":"syn4f"}).sort_index().fillna(0))