from collections import defaultdict
from typing import Set, Union, Dict, Iterable

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA

from pymutspec.constants import *


class CodonAnnotation:
    nucl_order = possible_nucls

    def __init__(self, gencode: Union[NCBICodonTableDNA, int]):
        self.codontable = self._prepare_codontable(gencode)
        self._syn_codons, self._ff_codons = self.__extract_syn_codons()
        self.possible_ff_contexts = self.__extract_possible_ff_contexts()
        self.possible_syn_contexts = self.__extract_possible_syn_contexts()
        self.startcodons, self.stopcodons = self.read_start_stop_codons(gencode)

    def is_fourfold(self, cdn: str):
        """Check if codon is neutral in 3rd position"""
        return cdn in self._ff_codons

    def translate_codon(self, cdn: str) -> str:
        """Translate codon to animo acid"""
        if isinstance(cdn, str):
            return self.codontable.forward_table.get(cdn, "*")
        else:
            return cdn

    def is_syn_mut(self, cdn1: str, cdn2: str):
        """
        Check if mutation in the codon is synonymous

        Arguments
        ---------
        cdn1: str
            Mutating codon (reference)
        cdn2: str
            Codon with mutation (alternative)

        Return
        ---------
            True if mutation is synonymous else False
        """
        if not isinstance(cdn1, str) or not isinstance(cdn2, str):
            return False
        return self.translate_codon(cdn1) == self.translate_codon(cdn2)

    def get_syn_codons(self, cdn: str, pic: int) -> Set[str]:
        """
        Get all possible synonymous codons

        Arguments
        ---------
        cdn: str
            Codon
        pic: int, [0, 1, 2]
            Position In Codon in that mutation occured

        Return
        -------
        syn_codons: Set[str]
            possible synonymous codons
        """
        assert 0 <= pic <= 2, "pic must be 0-based and less than 3"
        syn_codons = self._syn_codons.get((cdn, pic), dict())
        return syn_codons

    def get_mut_type(self, cdn1: str, cdn2: str, pic: int):
        """
        Arguments
        ---------
        cdn1: str
            Mutating codon (reference)
        cdn2: str
            Codon with mutation (alternative)
        pic: int, [0, 1, 2]
            Position In Codon in that mutation occured

        Return
        ------
        (label, aa1, aa2)

        aa - amino acid

        label variants:
        - -3 - stopcodon to stopcodon
        - -2 - stopcodon loss
        - -1 - stopcodon gain
        -  0 - non synonymous sbs
        -  1 - synonymous sbs
        -  2 - synonymous fourfold sbs
        """
        if not isinstance(cdn1, str) or not isinstance(cdn2, str):
            return 0, None, None
            
        _where_mut = [int(x != y) for x, y in zip(cdn1, cdn2)]
        if pic < 0 or pic > 2:
            raise ValueError("Position in codon (pic) must be 0-based and less than 3")
        if sum(_where_mut) != 1 or _where_mut[pic] != 1:
            raise ValueError("One mutation must be in codons and it must be on position 'pic'")
        aa1 = self.translate_codon(cdn1)
        aa2 = self.translate_codon(cdn2)
        if aa1 == "*" and aa2 == "*":
            label = -3
        elif aa1 == "*" and aa2 != "*":
            label = -2
        elif aa1 != "*" and aa2 == "*":
            label = -1
        elif aa1 == aa2:
            label = 1
            if pic == 2 and self.is_fourfold(cdn1):
                label = 2
        else:
            label = 0

        return label, aa1, aa2

    def extract_mutations_simple(self, g1: np.ndarray, g2: np.ndarray):
        """
        Extract alterations of genome2 (g2) comparing to genome1 (g1)

        Conditions of mutation collecting:
        - only sbs must be in one codon
        - only sbs must be in one trinucleotide context of sbs
        - pentanucleotide context must contain only explicit nucleotides

        Arguments
        ---------
        g1: Iterable
            reference sequence (parent node)
        g2: Iterable
            alternative sequence (child node)

        Return
        ---------
        mut_df: pd.DataFrame
            table of mutations with columns:
            - Mut
            - Label
            - PosInGene
            - PosInCodon
            - RefCodon
            - AltCodon
            - RefAa
            - AltAa
        """
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        mutations = []
        # pass initial codon and last nucleotide without right context
        for pos in range(3, n - 1):
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
                continue  # pentanucleotide context must contain only explicit nucleotides

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

    def collect_exp_mut_freqs(self, cds: Union[str, Iterable[str]], mask: Iterable[Union[int, bool]] = None, labels=["all", "syn", "ff"]):
        """
        Calculate potential expected mutation counts for nucleotides and trinucleotides (context) 
        in cds gene

        Arguments
        ---------
        cds: string or iterable of strings, if len is not divisible by 3, last codon is not used in syn, syn4f and pos3 modes
            cds sequence with codon structure; 
        labels: List of label strings
            label could be one of ["all", "syn", "ff", "pos3"];
        mask:
            iterable that mask invariant positions in the cds;

        Return
        ---------
        sbs12_freqs: Dict[label, Dict[nucl, count]]
            for each label collected expected single nucleotide substitutions frequencies without contexts
        sbs192_freqs: Dict[label, Dict[context, count]]
            for each label collected expected single nucleotide substitutions frequencies with contexts
        """
        n = len(cds)
        if mask is not None and len(mask) != n:
            raise ValueError("Mask must have same lenght as cds")

        labels = set(labels)
        sbs12_freqs = {lbl: defaultdict(int) for lbl in labels}
        sbs192_freqs = {lbl: defaultdict(int) for lbl in labels}

        for pos in range(1, n - 1):
            if mask is not None and not mask[pos]:
                continue
            pic = pos % 3
            nuc = cds[pos]
            cdn = cds[pos - pic: pos - pic + 3]
            cdn = cdn if isinstance(cdn, str) else "".join(cdn)
            cxt = cds[pos - 1: pos + 2]
            cxt = cxt if isinstance(cxt, str) else "".join(cxt)
            sbs12_pattern = nuc + ">" + "{}"
            sbs192_pattern = cxt[0] + "[" + nuc + ">{}]" + cxt[-1]

            if ("syn" in labels or "syn_c" in labels) and len(cdn) == 3:
                syn_codons = self.get_syn_codons(cdn, pic)
                if "syn" in labels:
                    for alt_cdn in syn_codons:
                        alt_nuc = alt_cdn[pic]
                        sbs12_freqs["syn"][sbs12_pattern.format(alt_nuc)] += 1
                        sbs192_freqs["syn"][sbs192_pattern.format(alt_nuc)] += 1
                if "syn_c" in labels and len(syn_codons) > 0:
                    for alt_nuc in self.nucl_order:
                        sbs12_freqs["syn_c"][sbs12_pattern.format(alt_nuc)] += 1
                        sbs192_freqs["syn_c"][sbs192_pattern.format(alt_nuc)] += 1
            
            for alt_nuc in self.nucl_order:
                if alt_nuc == nuc:
                    continue
                if "all" in labels:
                    sbs12_freqs["all"][sbs12_pattern.format(alt_nuc)] += 1
                    sbs192_freqs["all"][sbs192_pattern.format(alt_nuc)] += 1
                if "pos3" in labels and pic == 2:
                    sbs12_freqs["pos3"][sbs12_pattern.format(alt_nuc)] += 1
                    sbs192_freqs["pos3"][sbs192_pattern.format(alt_nuc)] += 1
                if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                    sbs12_freqs["ff"][sbs12_pattern.format(alt_nuc)] += 1
                    sbs192_freqs["ff"][sbs192_pattern.format(alt_nuc)] += 1

        return sbs12_freqs, sbs192_freqs

    def collect_exp_muts(self, cds, mask=None, labels=["syn"]):
        """
        Calculate potential expected mutations in cds gene

        Arguments
        ---------
        cds: string or iterable of strings, if len is not divisible by 3, last codon is not used in syn, syn4f and pos3 modes
            cds sequence with codon structure; 
        labels: List of label strings
            label could be one of ["all", "syn", "ff", "pos3"];
        mask:
            iterable that mask invariant positions in the cds;

        Return
        ---------
            sbs_table: pd.DataFrame
        """
        n = len(cds)
        if mask is not None and len(mask) != n:
            raise ValueError("Mask must have same lenght as cds")

        labels = set(labels)
        data = []
        for pos in range(1, n - 1):
            if mask is not None and not mask[pos]:
                continue
            pic = pos % 3
            nuc = cds[pos]
            cdn = cds[pos - pic: pos - pic + 3]
            cdn = cdn if isinstance(cdn, str) else "".join(cdn)
            cxt = cds[pos - 1: pos + 2]
            cxt = cxt if isinstance(cxt, str) else "".join(cxt)
            sbs192_pattern = cxt[0] + "[" + nuc + ">{}]" + cxt[-1]

            if ("syn" in labels or "syn_c" in labels) and len(cdn) == 3:
                syn_codons = self.get_syn_codons(cdn, pic)
                if "syn" in labels:
                    for alt_cdn in syn_codons:
                        alt_nuc = alt_cdn[pic]
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "syn",
                        })
                if "syn_c" in labels and len(syn_codons) > 0:
                    for alt_nuc in self.nucl_order:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "syn_c",
                        })

            for alt_nuc in self.nucl_order:
                if alt_nuc == nuc:
                    continue
                if "all" in labels:
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": sbs192_pattern.format(alt_nuc),
                        "Cdn": cdn, "Label": "all",
                    })
                if "pos3" in labels and pic == 2:
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": sbs192_pattern.format(alt_nuc),
                        "Cdn": cdn, "Label": "pos3",
                    })
                if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": sbs192_pattern.format(alt_nuc),
                        "Cdn": cdn, "Label": "syn4f",
                    })

        exp_sbs = pd.DataFrame(data)
        return exp_sbs

    def __extract_syn_codons(self):
        """
        extract synonymous (codons that mutate without amino acid change) 
        and fourfold codons from codon table  TODO rewrite

        Used in expected mutspec calculation

        Return
        -------
        syn_codons: Dict[Tuple(str, int), Set[str]]
            mapping (codon, pos_in_codon) --> set of possible syn codons
        ff_codons: set[str]
            set of ff codons (neutral on 3rd position)
        """
        aa2codons = defaultdict(set)
        for cdn, aa in self.codontable.forward_table.items():
            aa2codons[aa].add(cdn)

        syn_codons = defaultdict(set)
        for aa, codons in aa2codons.items():
            if len(codons) > 1:
                interim_dct = defaultdict(set)
                for i, slc in enumerate([slice(1, 3), slice(0, 3, 2), slice(0, 2)]):
                    for cdn in codons:
                        cdn_stump = cdn[slc]
                        interim_dct[(cdn_stump, i)].add(cdn)

                for key, aa_syn_codons in interim_dct.items():
                    if len(aa_syn_codons) > 1:
                        pic = key[1]
                        for cdn1 in aa_syn_codons:
                            for cdn2 in aa_syn_codons:
                                if cdn1 != cdn2:
                                    syn_codons[(cdn1, pic)].add(cdn2)
        ff_codons = set()
        for (cdn, pic), codons in syn_codons.items():
            if len(codons) == 3 and pic == 2:
                ff_codons.add(cdn)
        return dict(syn_codons), ff_codons

    def __extract_possible_ff_contexts(self) -> Set[str]:
        "Extract all contexts of neutral fourfold positions for current genetic code"
        possible_ff_contexts = set()
        for cdn in self._ff_codons:
            stump = cdn[1:]
            for nucl in self.nucl_order:
                cxt = stump + nucl
                possible_ff_contexts.add(cxt)
        return possible_ff_contexts

    def __extract_possible_syn_contexts(self) -> Set[str]:
        "Extract all contexts of neutral fourfold positions for current genetic code"
        possible_syn_contexts = set()
        for cdn, pic in self._syn_codons:
            if pic == 1:
                possible_syn_contexts.add(cdn)
                continue
            elif pic == 0:
                stump = cdn[:-1]
                for nucl in self.nucl_order:
                    cxt = nucl + stump
                    possible_syn_contexts.add(cxt)
            elif pic == 2:
                stump = cdn[1:]
                for nucl in self.nucl_order:
                    cxt = stump + nucl
                    possible_syn_contexts.add(cxt)
            else:
                raise RuntimeError()
        return possible_syn_contexts

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
                    cdn = [self.nucl_order[x] for x in [i, j, k]]
                    yield cdn, codon_proba

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
                    cdn = tuple(self.nucl_order[_] for _ in (i, j, k))

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
                            yield cdn, mut_context, full_proba
                    else:
                        yield cdn, cdn, codon_proba

    @staticmethod
    def read_start_stop_codons(codontable: Union[NCBICodonTableDNA, int]):
        codontable = CodonAnnotation._prepare_codontable(codontable)
        return set(codontable.start_codons), set(codontable.stop_codons)


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

    Return
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

    if "syn" in pivot_mutations.columns and "syn4f" in pivot_mutations.columns:
        pivot_mutations["syn"] += pivot_mutations["syn4f"]
    return pivot_mutations
