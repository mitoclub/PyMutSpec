from collections import defaultdict
from typing import Set, Union, Dict

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA

from .constants import *


class CodonAnnotation:
    nucl_order = possible_nucls

    def __init__(self, gencode: Union[NCBICodonTableDNA, int]):
        self.codontable = self._prepare_codontable(gencode)
        self._syn_codons, self._ff_codons = self.__extract_syn_codons()
        self.possible_ff_contexts = self.__extract_possible_ff_contexts()
        self.startcodons, self.stopcodons = self.read_start_stop_codons(gencode)

    def is_four_fold(self, codon):
        return codon in self._ff_codons

    def get_aa(self, codon: str):
        return self.codontable.forward_table.get(codon, "*")

    def is_syn_codons(self, codon1: str, codon2: str):
        if not isinstance(codon1, str) or not isinstance(codon2, str):
            return False
        return self.get_aa(codon1) == self.get_aa(codon2)

    def get_syn_number(self, cdn: str, pic: int):
        """return number of possible syn codons"""
        assert 0 <= pic <= 2, "pic must be 0-based and less than 3"
        return self._syn_codons.get((cdn, pic), 0)

    def get_mut_type(self, codon1: str, codon2: str, pic: int):
        """
        returned label variants:
        - -3 - stopcodon to stopcodon
        - -2 - stopcodon loss
        - -1 - stopcodon gain
        -  0 - non synonymous sbs
        -  1 - synonymous sbs
        -  2 - synonymous fourfold sbs

        return (label, aa1, aa2)
        """
        assert codon1 != codon2, "codons must be different"
        assert 0 <= pic <= 2, "pic must be 0-based and less than 3"
        aa1 = self.get_aa(codon1)
        aa2 = self.get_aa(codon2)
        if aa1 == "*" and aa2 == "*":
            label = -3
        elif aa1 == "*" and aa2 != "*":
            label = -2
        elif aa1 != "*" and aa2 == "*":
            label = -1
        elif aa1 == aa2:
            label = 1
            if pic == 2 and self.is_four_fold(codon1):
                label = 2
        else:
            label = 0

        return label, aa1, aa2

    def extract_mutations_simple(self, g1: np.ndarray, g2: np.ndarray):
        """
        Extract alterations of g2 comparing to g1

        params:
        - g1 - reference sequence (parent node)
        - g2 - alternative sequence (child node)

        conditions:
        - in one codon could be only sbs
        - in the context of one mutation couldn't be other sbs
        - indels are not sbs and codons and contexts with sbs are not considered

        return:
        - dataframe of mutations
        """
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        mutations = []
        for pos in range(1, n - 1):
            pic = pos % 3  # 0-based
            cdn_start = pos - pic
            cdn1 = g1[cdn_start: cdn_start + 3]
            cdn2 = g2[cdn_start: cdn_start + 3]
            mut_cxt1 = g1[pos - 1: pos + 2]
            mut_cxt2 = g2[pos - 1: pos + 2]
            cdn1_str = "".join(cdn1)
            cdn2_str = "".join(cdn2)

            up_nuc1, up_nuc2 = mut_cxt1[0], mut_cxt2[0]
            nuc1, nuc2 = mut_cxt1[1], mut_cxt2[1]
            down_nuc1, down_nuc2 = mut_cxt1[2], mut_cxt2[2]

            if nuc1 == nuc2 or up_nuc1 != up_nuc2 or down_nuc1 != down_nuc2:
                continue
            if sum([cdn1[_] == cdn2[_] for _ in range(3)]) != 2:
                continue
            if len(set(g1[pos - 2: pos + 3]).union(g2[pos - 2: pos + 3]) - set(self.nucl_order)) != 0:
                continue

            label, aa1, aa2 = self.get_mut_type(cdn1_str, cdn2_str, pic)
            sbs = {
                "Mut": f"{up_nuc1}[{nuc1}>{nuc2}]{down_nuc1}",
                "Label": np.int8(label),
                "PosInGene": np.int32(pos + 1),
                "PosInCodon": np.int8(pic + 1),
                "RefCodon": cdn1_str,
                "AltCodon": cdn2_str,
                "RefAa": aa1,
                "AltAa": aa2,
            }
            mutations.append(sbs)

        mut_df = pd.DataFrame(mutations)
        return mut_df

    def collect_obs_mut_freqs(self, genome: np.ndarray):
        n = len(genome)
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        nucl_freqs = {lbl: defaultdict(int) for lbl in ("all", "syn", "ff")}
        cxt_freqs = {lbl: defaultdict(int) for lbl in ("all", "syn", "ff")}

        for pos in range(1, n - 1):
            pic = pos % 3
            nuc = genome[pos]
            cdn = genome[pos - pic: pos - pic + 3]
            cxt = genome[pos - 1: pos + 2]
            cdn_str = "".join(cdn)
            cxt_str = "".join(cxt)

            nucl_freqs["all"][nuc] += 1
            cxt_freqs["all"][cxt_str] += 1

            syn_num = self.get_syn_number(cdn_str, pic)
            if syn_num > 0:
                nucl_freqs["syn"][nuc] += syn_num
                cxt_freqs["syn"][cxt_str] += syn_num

                if pic == 2 and self.is_four_fold(cdn_str):
                    nucl_freqs["ff"][nuc] += syn_num
                    cxt_freqs["ff"][cxt_str] += syn_num

        return nucl_freqs, cxt_freqs

    def __extract_syn_codons(self):
        """ extract synonymous (codons that mutate without amino acid change) 
        and fourfold codons from codon table

        usefull function for expected mutspec (filtration)

        return mapping[(cdn, pic)] of syn codons and set of ff codons
        """
        aa2codons = defaultdict(set)
        for codon, aa in self.codontable.forward_table.items():
            aa2codons[aa].add(codon)

        syn_codons = defaultdict(int)
        for aa, codons in aa2codons.items():
            if len(codons) > 1:
                interim_dct = defaultdict(set)
                for i, slc in enumerate([slice(1, 3), slice(0, 3, 2), slice(0, 2)]):
                    for codon in codons:
                        cdn_stump = codon[slc]
                        interim_dct[(cdn_stump, i)].add(codon)

                for key, aa_syn_codons in interim_dct.items():
                    if len(aa_syn_codons) > 1:
                        pic = key[1]
                        for cdn in aa_syn_codons:
                            syn_codons[(cdn, pic)] += len(aa_syn_codons) - 1
        ff_codons = set()
        for (cdn, pic), num in syn_codons.items():
            if num == 3 and pic == 2:
                ff_codons.add(cdn)
        return dict(syn_codons), ff_codons

    def __extract_possible_ff_contexts(self) -> Set[str]:
        possible_ff_contexts = set()
        for cdn in self._ff_codons:
            stump = cdn[1:]
            for nucl in self.nucl_order:
                cxt = stump + nucl
                possible_ff_contexts.add(cxt)
        return possible_ff_contexts

    @staticmethod
    def _prepare_codontable(codontable: Union[NCBICodonTableDNA, int]):
        if isinstance(codontable, NCBICodonTableDNA):
            pass
        elif isinstance(codontable, int):
            codontable = CodonTable.unambiguous_dna_by_id[codontable]
        else:
            ValueError("passed codontable is not appropriate")
        return codontable

    def _codon_iterator(self, codon_states: np.ndarray, cutoff=0.01):
        assert codon_states.shape == (3, 4)
        for i, p1 in enumerate(codon_states[0]):
            if p1 < cutoff:
                continue
            for j, p2 in enumerate(codon_states[1]):
                if p2 < cutoff:
                    continue
                for k, p3 in enumerate(codon_states[2]):
                    if p3 < cutoff:
                        continue
                    codon_proba = p1 * p2 * p3
                    if codon_proba < cutoff:
                        continue
                    codon = [self.nucl_order[x] for x in [i, j, k]]
                    yield codon, codon_proba

    def _sample_context(self, pos, pic, genome: np.ndarray, cutoff=0.01):
        nuc_cutoff = cutoff * 5
        codon_states = genome[pos - pic: pos - pic + 3]
        extra_codon_states = genome[pos + pic - 1]  # doesn't mean if pic == 1
        # gaps are not appropriate
        if np.any(codon_states.sum(1) == 0) or extra_codon_states.sum() == 0:
            return

        for i, p1 in enumerate(codon_states[0]):
            if p1 < nuc_cutoff:
                continue
            for j, p2 in enumerate(codon_states[1]):
                if p2 < nuc_cutoff:
                    continue
                for k, p3 in enumerate(codon_states[2]):
                    if p3 < nuc_cutoff:
                        continue
                    codon_proba = p1 * p2 * p3
                    if codon_proba < cutoff:
                        continue
                    codon = tuple(self.nucl_order[_] for _ in (i, j, k))

                    if pic != 1:
                        for m, p4 in enumerate(extra_codon_states):
                            if p4 < nuc_cutoff:
                                continue
                            full_proba = codon_proba * p4
                            if full_proba < cutoff:
                                continue
                            if pic == 0:
                                mut_context = tuple(self.nucl_order[_] for _ in (m, i, j))
                            elif pic == 2:
                                mut_context = tuple(self.nucl_order[_] for _ in (j, k, m))
                            yield codon, mut_context, full_proba
                    else:
                        yield codon, codon, codon_proba

    @staticmethod
    def read_start_stop_codons(codontable: Union[NCBICodonTableDNA, int]):
        codontable = CodonAnnotation._prepare_codontable(codontable)
        return set(codontable.start_codons), set(codontable.stop_codons)


def calculate_mutspec(mutations: pd.DataFrame, freqs: Dict[str, float], label: str, use_context=True, use_proba=True):
    """
    mutations dataframe must contain 2 columns:
    - Mut (X[N1>N2]Y)
    - Label (-1, 0, 1, 2)
    - ProbaFull (optional, only for use_proba=True)
    """
    mut = mutations.copy()
    _cols = ["Label", "Mut", "ProbaFull"] if use_proba else ["Label", "Mut"]
    for c in _cols:
        assert c in mut.columns, f"Column {c} is not in mut df"

    labels = {"syn", "ff", "all"}
    if isinstance(label, str):
        label = label.lower()
        if label not in labels:
            raise ValueError(f"pass the appropriate label: {labels}")
        if label == "syn":
            label = 1
        elif label == "ff":
            label = 2
        elif label == "all":
            label = 0
    else:
        raise ValueError(f"pass the appropriate label: {labels}")

    if use_context:
        col_mut = "Mut"
        full_sbs = possible_sbs192_set
    else:
        mut["MutBase"] = mut["Mut"].str.slice(2, 5)
        col_mut = "MutBase"
        full_sbs = possible_sbs12_set

    if use_proba:
        mutspec = mut[mut["Label"] >= label].groupby(col_mut)["ProbaFull"].sum().reset_index()
    else:
        mutspec = mut[mut["Label"] >= label][col_mut].value_counts().reset_index()
    mutspec.columns = ["Mut", "ObsFr"]

    # fill unobserved mutations by zeros
    mutspec_appendix = []
    unobserved_sbs = full_sbs.difference(mutspec["Mut"].values)
    for usbs in unobserved_sbs:
        mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
    mutspec = pd.concat([mutspec, pd.DataFrame(mutspec_appendix)], ignore_index=True)

    if use_context:
        sbs = mutspec["Mut"]
        mutspec["Context"] = sbs.str.get(0) + sbs.str.get(2) + sbs.str.get(-1)
    else:
        mutspec["Context"] = mutspec["Mut"].str.get(0)

    mutspec["ExpFr"] = mutspec["Context"].map(freqs)
    mutspec["RawMutSpec"] = (mutspec["ObsFr"] / mutspec["ExpFr"]).fillna(0)
    mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
    mutspec.drop("Context", axis=1, inplace=True)
    return mutspec


def mutations_summary(mutations: pd.DataFrame, gene_col=None, proba_col=None, gene_name_mapper: dict = None):
    """
    form mutations annotation: how many synonymous, fourfold or stop loss/gain observed in the table

    Arguments
    ---------
    mutations: pd.DataFrame 
        table must contain at least 2 columns:
        - Mut: str; Pattern: '[ACGT]\[[ACGT]>[ACGT]\][ACGT]'
        - Label: int; [-3, 2]. See CodonAnnotation.get_mut_type
        - $gene_col, optional. If gene_col=None annotation will be formed on full mutations without genes splitting
        - $proba_col, optional. If proba_col=None each row of table assumed to be mutation, else probabilities will be used

    gene_col: str
        Column containing gene name in that mutation was observed
    proba_col: str
        Column containing probability of mutations
    gene_name_mapper: dict
        mapping for gene names. Use when gene_col contains indexses of genes

    Returns
    -------
    pivot_mutations: pd.DataFrame 
        table with mutations annotation
    """
    mutations = mutations.copy()
    mut_pattern = "[ACGT]\[[ACGT]>[ACGT]\][ACGT]"
    label_mapper = {
        -3: "6Stop to stop",
        -2: "4Stop loss",
        -1: "5Stop gain",
        0: "1non-syn",
        1: "2syn",
        2: "3syn4f",
    }
    grp = ["Label"]
    if gene_col is not None:
        grp.append(gene_col)

    if proba_col is None:
        proba_col = "_ProbaFull"
        mutations[proba_col] = 1
    assert proba_col in mutations.columns

    mutations_descr = mutations[
        (mutations.Mut.str.fullmatch(mut_pattern))
    ].groupby(grp)[proba_col].sum().reset_index()

    mutations_descr["Label"] = mutations_descr.Label.map(label_mapper)
    pivot_mutations = mutations_descr.pivot_table(proba_col, gene_col, "Label", fill_value=0)
    pivot_mutations.columns = [x[1:] for x in pivot_mutations.columns]

    if gene_name_mapper is not None:
        pivot_mutations.index = pivot_mutations.index.map(gene_name_mapper)
    return pivot_mutations


translator = str.maketrans("ACGT", "TGCA")

def rev_comp(mut: str):
    new_mut = mut[-1] + mut[1:-1] + mut[0]
    new_mut = new_mut.translate(translator)
    return new_mut
