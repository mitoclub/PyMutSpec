import os
import sys
from collections import defaultdict
from typing import Set, Union, Dict, Iterable
import multiprocessing as mp
from warnings import warn

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA
from ete3 import PhyloTree

from pymutspec.constants import *
from pymutspec.utils import basic_logger
from ..io import GenomeStatesTotal
from .tree import iter_tree_edges, calc_phylocoefs

nucl_order5 = ['A', 'C', 'G', 'T', '-']


class CodonAnnotation:
    nucl_order = possible_nucls

    def __init__(self, gencode: Union[NCBICodonTableDNA, int]):
        self.codontable = self._prepare_codontable(gencode)
        self._syn_codons, self._ff_codons = self.__extract_syn_codons()
        self.possible_ff_contexts  = self.__extract_possible_ff_contexts()
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
        nucleotides = set("ACGTacgt")

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
            if len(set(g1[pos -1: pos + 2]).union(g2[pos - 1: pos + 2]) - nucleotides) != 0:
                continue  # trinucleotide context must contain only explicit nucleotides

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
        
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

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
            syn_codons = self.get_syn_codons(cdn, pic)

            if "syn" in labels:
                for alt_cdn in syn_codons:
                    alt_nuc = alt_cdn[pic]
                    sbs12_freqs["syn"][sbs12_pattern.format(alt_nuc)] += 1
                    sbs192_freqs["syn"][sbs192_pattern.format(alt_nuc)] += 1
            if "nonsyn" in labels:
                syn_alt_nucs = [cdn[pic] for cdn in syn_codons]
                syn_alt_nucs.append(nuc)
                nonsyn_alt_nucs = set(self.nucl_order).difference(syn_alt_nucs)
                for alt_nuc in nonsyn_alt_nucs:
                    sbs12_freqs["nonsyn"][sbs12_pattern.format(alt_nuc)] += 1
                    sbs192_freqs["nonsyn"][sbs192_pattern.format(alt_nuc)] += 1

            for alt_nuc in self.nucl_order:
                if alt_nuc == nuc:
                    continue
                cur_sbs12  = sbs12_pattern.format(alt_nuc)
                cur_sbs192 = sbs192_pattern.format(alt_nuc)

                if "all" in labels:
                    sbs12_freqs["all"][cur_sbs12] += 1
                    sbs192_freqs["all"][cur_sbs192] += 1
                if "pos3" in labels and pic == 2:
                    sbs12_freqs["pos3"][cur_sbs12] += 1
                    sbs192_freqs["pos3"][cur_sbs192] += 1
                if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                    sbs12_freqs["ff"][cur_sbs12] += 1
                    sbs192_freqs["ff"][cur_sbs192] += 1
                if "syn_c" in labels and len(syn_codons) > 0:
                    sbs12_freqs["syn_c"][cur_sbs12] += 1
                    sbs192_freqs["syn_c"][cur_sbs192] += 1

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

        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

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
            syn_codons = self.get_syn_codons(cdn, pic)

            if "syn" in labels:
                for alt_cdn in syn_codons:
                    alt_nuc = alt_cdn[pic]
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": sbs192_pattern.format(alt_nuc),
                        "Cdn": cdn, "Label": "syn",
                    })
            if "nonsyn" in labels:
                syn_alt_nucs = [cdn[pic] for cdn in syn_codons]
                syn_alt_nucs.append(nuc)
                nonsyn_alt_nucs = set(self.nucl_order).difference(syn_alt_nucs)
                for alt_nuc in nonsyn_alt_nucs:
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": sbs192_pattern.format(alt_nuc),
                        "Cdn": cdn, "Label": "nonsyn",
                    })

            for alt_nuc in self.nucl_order:
                if alt_nuc == nuc:
                    continue
                cur_sbs192 = sbs192_pattern.format(alt_nuc)

                if "all" in labels:
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": cur_sbs192,
                        "Cdn": cdn, "Label": "all",
                    })
                if "pos3" in labels and pic == 2:
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": cur_sbs192,
                        "Cdn": cdn, "Label": "pos3",
                    })
                if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": cur_sbs192,
                        "Cdn": cdn, "Label": "syn4f",
                    })
                if "syn_c" in labels and len(syn_codons) > 0:
                    data.append({
                        "Pos": pos + 1, "Pic": pic + 1,
                        "Mut": cur_sbs192,
                        "Cdn": cdn, "Label": "syn_c",
                    })

        exp_sbs = pd.DataFrame(data)
        return exp_sbs

    def extract_mutations_proba(self, g1: np.ndarray, g2: np.ndarray, phylocoef: float, mut_proba_cutoff: float):
        """
        Extract alterations of g2 comparing to g1
        TODO

        conditions:
        - in one codon could be only sbs
        - in the context of one mutation couldn't be other sbs
        - indels are not sbs and codons and contexts with sbs are not considered

        Arguments
        - g1 - reference sequence (parent node)
        - g2 - alternative sequence (child node)

        return:
        - mut - dataframe of mutations
        """
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        mutations = []
        # pass initial codon and last nucleotide without right context
        for pos in range(3, n - 1):
            pic = pos % 3  # 0-based
            for cdn1, mut_cxt1, proba1 in self.sample_context(pos, pic, g1, mut_proba_cutoff / phylocoef):
                cdn1_str = "".join(cdn1)
                for cdn2, mut_cxt2, proba2 in self.sample_context(pos, pic, g2, mut_proba_cutoff / phylocoef):
                    p_adj = proba1 * proba2 * phylocoef
                    if p_adj < mut_proba_cutoff:
                        continue

                    cdn2_str = "".join(cdn2)

                    up_nuc1, up_nuc2 = mut_cxt1[0], mut_cxt2[0]
                    nuc1, nuc2 = mut_cxt1[1], mut_cxt2[1]
                    down_nuc1, down_nuc2 = mut_cxt1[2], mut_cxt2[2]

                    if nuc1 == nuc2 or up_nuc1 != up_nuc2 or down_nuc1 != down_nuc2:
                        continue
                    if sum([cdn1[_] == cdn2[_] for _ in range(3)]) != 2:
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
                        "ProbaRef": proba1,
                        "ProbaMut": proba1 * proba2,
                        "ProbaFull": p_adj,
                    }
                    mutations.append(sbs)

        mut_df = pd.DataFrame(mutations)
        return mut_df

    def collect_exp_mut_freqs_proba(
            self, cds: np.ndarray, phylocoef: float, 
            mask: Iterable[Union[int, bool]] = None, 
            labels = ["all", "syn", "ff"],  mut_proba_cutoff=0.05,
        ):
        n = len(cds)
        if mask is not None and len(mask) != n:
            msg = f"Mask (len = {len(mask)}) must have same lenght as cds (len = {n})"
            print(msg, file=sys.stderr)
            # logger.error(msg)
            # logger.info("Termination")
            raise ValueError(msg)

        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        assert 0 < phylocoef <= 1, "Evol coefficient must be between 0 and 1"

        labels = set(labels)
        sbs12_freqs = {lbl: defaultdict(int) for lbl in labels}
        sbs192_freqs = {lbl: defaultdict(int) for lbl in labels}

        for pos in range(1, n - 1):
            if mask is not None and not mask[pos]:
                continue
            pic = pos % 3  # 0-based
            for cdn_tuple, cxt, p in self.sample_context(pos, pic, cds, mut_proba_cutoff / phylocoef):
                # we don't use low-probability mutations by unite cutoff `mut_proba_cutoff / phylocoef`
                cdn = "".join(cdn_tuple)
                p_adj = p * phylocoef  # adjusted probability
                nuc = cxt[1]
                sbs12_pattern = nuc + ">" + "{}"
                sbs192_pattern = cxt[0] + "[" + nuc + ">{}]" + cxt[-1]
                syn_codons = self.get_syn_codons(cdn, pic)

                if "syn" in labels:
                    for alt_cdn in syn_codons:
                        alt_nuc = alt_cdn[pic]
                        sbs12_freqs["syn"][sbs12_pattern.format(alt_nuc)] += p_adj
                        sbs192_freqs["syn"][sbs192_pattern.format(alt_nuc)] += p_adj
                if "nonsyn" in labels:
                    syn_alt_nucs = [cdn[pic] for cdn in syn_codons]
                    syn_alt_nucs.append(nuc)
                    nonsyn_alt_nucs = set(self.nucl_order).difference(syn_alt_nucs)
                    for alt_nuc in nonsyn_alt_nucs:
                        sbs12_freqs["nonsyn"][sbs12_pattern.format(alt_nuc)] += p_adj
                        sbs192_freqs["nonsyn"][sbs192_pattern.format(alt_nuc)] += p_adj

                for alt_nuc in self.nucl_order:
                    if alt_nuc == nuc:
                        continue
                    cur_sbs12  = sbs12_pattern.format(alt_nuc)
                    cur_sbs192 = sbs192_pattern.format(alt_nuc)

                    if "all" in labels:
                        sbs12_freqs["all"][cur_sbs12] += p_adj
                        sbs192_freqs["all"][cur_sbs192] += p_adj
                    if "pos3" in labels and pic == 2:
                        sbs12_freqs["pos3"][cur_sbs12] += p_adj
                        sbs192_freqs["pos3"][cur_sbs192] += p_adj
                    if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                        sbs12_freqs["ff"][cur_sbs12] += p_adj
                        sbs192_freqs["ff"][cur_sbs192] += p_adj
                    if "syn_c" in labels and len(syn_codons) > 0:
                        sbs12_freqs["syn_c"][cur_sbs12] += p_adj
                        sbs192_freqs["syn_c"][cur_sbs192] += p_adj
                        
        return sbs12_freqs, sbs192_freqs

    def collect_exp_muts_proba(
            self, cds: np.ndarray, phylocoef: float, 
            mask: Iterable[Union[int, bool]] = None, 
            labels = ["all", "syn", "ff"],  mut_proba_cutoff=0.05,
        ):
        n = len(cds)
        if mask is not None and len(mask) != n:
            msg = f"Mask (len = {len(mask)}) must have same lenght as cds (len = {n})"
            print(msg, file=sys.stderr)
            # logger.error(msg)
            # logger.info("Termination")
            raise ValueError(msg)

        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        assert 0 < phylocoef <= 1, "Evol coefficient must be between 0 and 1"

        labels = set(labels)
        data = []
        for pos in range(1, n - 1):
            if mask is not None and not mask[pos]:
                continue
            pic = pos % 3  # 0-based
            for cdn_tuple, cxt, p in self.sample_context(pos, pic, cds, mut_proba_cutoff / phylocoef):
                # we don't use low-probability mutations by unite cutoff `mut_proba_cutoff / phylocoef`
                cdn = "".join(cdn_tuple)
                p_adj = p * phylocoef  # adjusted probability
                nuc = cxt[1]
                sbs192_pattern = cxt[0] + "[" + nuc + ">{}]" + cxt[-1]
                syn_codons = self.get_syn_codons(cdn, pic)

                if "syn" in labels:
                    for alt_cdn in syn_codons:
                        alt_nuc = alt_cdn[pic]
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "syn",
                            "Proba": p_adj,
                        })
                if "nonsyn" in labels:
                    syn_alt_nucs = [cdn[pic] for cdn in syn_codons]
                    syn_alt_nucs.append(nuc)
                    nonsyn_alt_nucs = set(self.nucl_order).difference(syn_alt_nucs)
                    for alt_nuc in nonsyn_alt_nucs:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "nonsyn",
                            "Proba": p_adj,
                        })

                for alt_nuc in self.nucl_order:
                    if alt_nuc == nuc:
                        continue
                    cur_sbs192 = sbs192_pattern.format(alt_nuc)

                    if "all" in labels:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": cur_sbs192,
                            "Cdn": cdn, "Label": "all",
                            "Proba": p_adj,
                        })
                    if "pos3" in labels and pic == 2:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": cur_sbs192,
                            "Cdn": cdn, "Label": "pos3",
                            "Proba": p_adj,
                        })
                    if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": cur_sbs192,
                            "Cdn": cdn, "Label": "syn4f",
                            "Proba": p_adj,
                        })
                    if "syn_c" in labels and len(syn_codons) > 0:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "syn_c",
                            "Proba": p_adj,
                        })
        exp_sbs = pd.DataFrame(data)
        return exp_sbs

    def sample_context(self, pos, pic, genome: np.ndarray, cutoff=0.05):
        codon_states = genome[pos - pic: pos - pic + 3]        
        extra_codon_states = genome[pos + pic - 1]  # doesn't mean if pic == 1
        # gaps are not appropriate
        if np.any(codon_states.sum(1) == 0) or extra_codon_states.sum() == 0:
            return

        a, b, c = codon_states
        probas = a * b[:, None] * c[:, None, None]
        _ii = 0  # increment index if there are 4th context nucl
        if pic != 1:
            probas = probas * extra_codon_states[:, None, None, None]
            _ii = 1

        indexes = np.where(probas > cutoff)
        for idx in range(len(indexes[0])):
            i, j, k = indexes[2+_ii][idx], indexes[1+_ii][idx], indexes[0+_ii][idx]
            m = indexes[0][idx]

            cdn = tuple(self.nucl_order[_] for _ in (i, j, k))
            if pic == 0:
                mut_context = tuple(self.nucl_order[_] for _ in (m, i, j))
                full_proba = probas[m, k, j, i]
            elif pic == 2:
                mut_context = tuple(self.nucl_order[_] for _ in (j, k, m))
                full_proba = probas[m, k, j, i]
            elif pic == 1:
                mut_context = cdn
                full_proba = probas[k, j, i]
            
            yield cdn, mut_context, full_proba

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


class Branch:
    def __init__(self, index: int, ref_name: str, alt_name: str, 
                 ref_genome: pd.DataFrame, alt_genome: pd.DataFrame, 
                 phylocoef: float, proba_cutoff: float, genome_size: int, 
                 site2index: Dict[int, int], logger=None) -> None:
        self.index = index
        self.ref_name = ref_name
        self.alt_name = alt_name
        self.ref_genome = ref_genome
        self.alt_genome = alt_genome
        self.phylocoef = phylocoef
        self.proba_cutoff = proba_cutoff
        self.genome_size = genome_size
        self.site2index = site2index
        self.logger = logger or basic_logger()

    def is_substitution(self, cxt1, cxt2):
        """
        check that contexts contain substitution.

        Example inputs and outputs:
        - (CAT,CGG) -> False (T != G)
        - (CAT,CAG) -> False (central nuclaotide are constant)
        - (ACTAC,ACAAC) -> True
        - (ACTAC,AC-AC) -> False (indel != substitution)
        """
        ok = True
        center_pos = len(cxt1) // 2
        if cxt1[center_pos] == cxt2[center_pos]:
            ok = False
        elif cxt1[center_pos] == '-' or cxt2[center_pos] == '-':
            ok = False
        elif cxt1[center_pos-1] != cxt2[center_pos-1] or \
             cxt1[center_pos+1] != cxt2[center_pos+1]:
            self.logger.warning(f"Neighbouring positions of substitution are different: {cxt1}, {cxt2}")
            ok = False
        return ok

    def process_branch(self):
        edge_mutations = self.extract_mutations_proba(
            self.ref_genome, self.alt_genome, self.phylocoef, self.proba_cutoff)
        
        edge_mutations["RefNode"] = self.ref_name
        edge_mutations["AltNode"] = self.alt_name

        mut_num = edge_mutations['ProbaFull'].sum() if 'ProbaFull' in edge_mutations.columns else len(edge_mutations)
        if mut_num == 0:
            self.logger.info(f"No mutations from branch {self.index:03} ({self.ref_name} - {self.alt_name})")
        elif mut_num > self.genome_size * 0.1:
            self.logger.warning(f"Observed too many mutations ({mut_num} > {self.genome_size} * 0.1) "
                           f"for branch ({self.ref_name} - {self.alt_name})")
        else:
            self.logger.info(f"{mut_num:.2f} mutations from branch {self.index:03} ({self.ref_name} - {self.alt_name})")

        return edge_mutations

    def extract_mutations_proba(self, g1: pd.DataFrame, g2: pd.DataFrame, 
                                phylocoef: float, mut_proba_cutoff: float):
        """
        Extract alterations of g2 comparing to g1

        conditions:
        - in one codon could be only sbs
        - in the context of one mutation couldn't be other sbs
        - indels are not sbs and codons and contexts with sbs are not considered

        Arguments
        - g1 - reference sequence (parent node)
        - g2 - alternative sequence (child node)

        return:
        - mut - dataframe of mutations
        """
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"

        mutations = []
        g1_most_prob = g1.values.argmax(1)
        g2_most_prob = g2.values.argmax(1)
        # pass initial and last nucleotides without context
        for site in g1.index[2:-2]:
            site_mutations = self.process_site(
                g1, g2, site, g1_most_prob, g2_most_prob, mut_proba_cutoff, phylocoef)
            mutations.extend(site_mutations)
        mut_df = pd.DataFrame(mutations)
        return mut_df

    def process_site(self, g1: pd.DataFrame, g2: pd.DataFrame, site: int,
                     g1_most_prob: np.ndarray, g2_most_prob: np.ndarray, 
                     mut_proba_cutoff: float, phylocoef: float):
        pos_in_states = self.site2index[site]
        states1, states2 = self.prepare_local_states(
            g1.values, g2.values, pos_in_states, g1_most_prob, g2_most_prob)
        if states1 is None:
            return []

        mutations = []
        for cxt1, proba1 in self.sample_context5(2, states1, mut_proba_cutoff):
            for cxt2, proba2 in self.sample_context5(2, states2, mut_proba_cutoff):
                if not self.is_substitution(cxt1, cxt2):
                    continue

                p = proba1 * proba2
                if p < mut_proba_cutoff:
                    continue
                
                p_adj = p * phylocoef
                
                cxt1 = ''.join(cxt1)
                cxt2 = ''.join(cxt2)
                muttype = f'{cxt1[0:2]}[{cxt1[2]}>{cxt2[2]}]{cxt1[3:]}'
                sbs = {
                    "Mut": muttype,
                    "Site": site,
                    "ProbaRef": proba1,
                    "ProbaMut": p,
                    "ProbaFull": p_adj,
                }
                mutations.append(sbs)
        return mutations

    def prepare_local_states_simple(self, g1: np.ndarray, g2: np.ndarray, pos: int):
        return g1[pos-2: pos+3], g2[pos-2: pos+3]

    def prepare_local_states(self, g1: np.ndarray, g2: np.ndarray, pos: int,
                             g1_most_prob: np.ndarray, g2_most_prob: np.ndarray):
        # prepare local states: exclude gaps from both seqs
        ## first, search upstream
        up_states1, up_states2 = [g1[pos]], [g2[pos]]
        down_states1, down_states2 = [], []
        for i in range(1, 10): # don't collect gaps at all
            cur_state1, cur_state2 = g1[pos-i], g2[pos-i]
            state_id1, state_id2 = g1_most_prob[pos-i], g2_most_prob[pos-i]
            if not (state_id1 == state_id2 == 4): # both states are not equal to gap
                up_states1.append(cur_state1)
                up_states2.append(cur_state2)
            if len(up_states1) == 3 or pos-i-1 <= 0:
                break

        ## second, search downstream
        for i in range(1, 10):
            cur_state1, cur_state2 = g1[pos+i], g2[pos+i]
            state_id1, state_id2 = g1_most_prob[pos+i], g2_most_prob[pos+i]
            if not (state_id1 == state_id2 == 4): # not equal to gap
                down_states1.append(cur_state1)
                down_states2.append(cur_state2)
            if len(down_states1) == 2 or pos+i+1 >= len(g1):
                break
        
        if any([len(up_states1) != 3, 
                len(up_states2) != 3, 
                len(down_states1) != 2 , 
                len(down_states2) != 2]):
            self.logger.warning(f'Cannot collect local contexts for position {pos}')
            return None, None

        up_states1.reverse()
        up_states2.reverse()
        up_states1.extend(down_states1)
        up_states2.extend(down_states2)
        states1 = np.array(up_states1)
        states2 = np.array(up_states2)

        return states1, states2

    def sample_context5(self, pos: int, g: np.ndarray, cutoff=0.25):
        probas = g[pos-2] * g[pos-1][:, None] * g[pos][:, None, None] * \
            g[pos+1][:, None, None, None] * g[pos+2][:, None, None, None, None]
        
        indexes = np.where(probas > cutoff)
        for idx in range(len(indexes[0])):
            cxt_ids = [indexes[i][idx] for i in list(range(len(indexes)))[::-1]]
            cxt = tuple(nucl_order5[x] for x in cxt_ids)
            pr = probas[tuple(cxt_ids[::-1])]
            yield cxt, pr

    def sample_context7(self, pos: int, g: np.ndarray, cutoff=0.25):
        probas = g[pos-3] * g[pos-2][:, None] * g[pos-1][:, None, None] * \
            g[pos][:, None, None, None] * g[pos+1][:, None, None, None, None] * \
                g[pos+2][:, None, None, None, None, None] * \
                    g[pos+3][:, None, None, None, None, None, None]
        
        indexes = np.where(probas > cutoff)
        for idx in range(len(indexes[0])):
            cxt_ids = [indexes[i][idx] for i in list(range(len(indexes)))[::-1]]
            cxt = tuple(nucl_order5[x] for x in cxt_ids)
            pr = probas[tuple(cxt_ids[::-1])]
            yield cxt, pr


class MutSpecExtractor(CodonAnnotation):
    def __init__(
            self, path_to_tree, outdir, gcode=2,
            use_proba=False, proba_cutoff=0.25, use_phylocoef=False,
            genomes: GenomeStatesTotal = None, num_processes=8, 
            logger=None,
        ):
        if not os.path.exists(path_to_tree):
            raise ValueError(f"Path to tree doesn't exist: '{path_to_tree}'")

        CodonAnnotation.__init__(self, gencode=gcode)

        self.num_processes = num_processes
        self.genomes = genomes
        self.gcode = gcode
        self.use_proba = use_proba
        self.proba_cutoff = proba_cutoff
        self.use_phylocoef = use_phylocoef if use_proba else False
        self.outdir = outdir
        self.logger = logger or basic_logger()
        
        self.logger.info(f"Using gencode {gcode}")
        self.logger.info(f"Use probabilities of genomic states: {use_proba}")
        self.logger.info(f"Use phylogenetic uncertainty coefficient: {use_phylocoef}")
        self.logger.info(f"Minimal probability for mutations to use: {proba_cutoff}")

        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        self.logger.info(
            f"Tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
        )
        rnd_genome = self.genomes.get_random_genome()
        self.logger.info(f"Number of MSA sites: {len(rnd_genome)}")

         # mapping of MSA sites to truncated states ids
        self.site2index = dict(zip(rnd_genome.index, range(len(rnd_genome))))

    def open_handles(self, outdir):
        self.handle = dict()
        self.handle["mut"]  = open(os.path.join(outdir, "mutations.tsv"), "w")
        self.logger.debug("Handles opened")

    def close_handles(self):
        for file in self.handle.values():
            file.close()
        self.logger.debug("Handles closed")

    def extract_mutspec_from_tree(self):
        t = self.num_processes
        if not isinstance(t, int) or t < 1:
            raise ValueError('num_processes must be positive integer')

        if t == 1:
            self._derive_mutspec()
        elif t > 1:
            self._derive_mutspec_parallel()

    def _derive_mutspec(self):
        self.logger.info("Start mutation extraction from tree")
        self.open_handles(self.outdir)
        add_header = defaultdict(lambda: True)
        total_mut_num = 0

        for edge_data in self.iter_branches():
            edge_mutations = self.process_branch(*edge_data[1:])
            mut_num = edge_mutations['ProbaFull'].sum() if self.use_proba and \
                'ProbaFull' in edge_mutations.columns else len(edge_mutations)
            total_mut_num += mut_num

            # dump current edge mutations
            self.dump_table(edge_mutations, self.handle["mut"], add_header["mut"])
            add_header["mut"] = False

        self.close_handles()
        self.logger.info(f"Processed {edge_data[1]} tree edges")
        self.logger.info(f"Observed {total_mut_num:.3f} substitutions")
        self.logger.info("Extraction of mutations from phylogenetic tree completed succesfully")

    def _derive_mutspec_parallel(self):
        self.logger.info("Start parallel mutation extraction from tree")

        with mp.Pool(processes=self.num_processes) as pool:
            genome_mutations_lst = pool.map(Branch.process_branch, self.iter_branches())

        genome_mutations = pd.concat(genome_mutations_lst)

        self.open_handles(self.outdir)
        self.dump_table(genome_mutations, self.handle["mut"], True)
        self.close_handles()

        total_mut_num = genome_mutations['ProbaFull'].sum() if self.use_proba and \
            'ProbaFull' in genome_mutations.columns else len(genome_mutations)
        self.logger.info(f"Observed {total_mut_num:.3f} substitutions")
        self.logger.info("Extraction of mutations from phylogenetic tree completed succesfully")

    def iter_branches(self):
        """yield (self, edge_id, ref_node_name, alt_node_name, ref_genome, alt_genome, phylocoef)"""
        # calculate phylogenetic uncertainty correction
        if self.use_phylocoef:
            phylocoefs = calc_phylocoefs(self.tree)

        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree), 1):
            if alt_node.name not in self.genomes.nodes:
                self.logger.warning(f"Skip edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{alt_node.name}' genome")
                continue
            if ref_node.name not in self.genomes.nodes:
                self.logger.warning(f"Skip edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{ref_node.name}' genome")
                continue

            # calculate phylogenetic uncertainty correction: 
            #   coef for ref edge * coef for alt edge
            if self.use_phylocoef:
                phylocoef = self.fp_format(phylocoefs[ref_node.name] * phylocoefs[alt_node.name])
            else:
                phylocoef = self.fp_format(1)

            # get genomes from storage
            ref_genome = self.genomes.get_genome(ref_node.name)
            alt_genome = self.genomes.get_genome(alt_node.name)
            self.logger.debug(f'Genomes retrieved for branch {ei}')

            yield Branch(ei, ref_node.name, alt_node.name, ref_genome, alt_genome, 
                         phylocoef, self.proba_cutoff, self.genomes.genome_size, self.site2index)

    def get_edges_data(self):
        return list(self.iter_branches())

    @staticmethod
    def dump_table(df: pd.DataFrame, handle, header=False):
        if header:
            handle.write("\t".join(df.columns) + "\n")
        handle.write(df.to_csv(sep="\t", index=None, header=None, float_format='%g'))

    def turn_to_MAP(self, states: np.ndarray):
        if isinstance(states, pd.DataFrame):
            warn('pd.DataFrame used as input instead of np.ndarray', ResourceWarning)
            states = states.values
        return ''.join([nucl_order5[x] for x in states.argmax(1)])


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
