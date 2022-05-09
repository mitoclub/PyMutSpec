from collections import defaultdict
from typing import Union, Dict

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA

from .constants import *


class CodonAnnotation:
    nucl_order = ["A", "C", "G", "T"]

    def __init__(self, codontable: Union[NCBICodonTableDNA, int], *args, **kwargs):
        self.codontable = self._prepare_codontable(codontable)
        self._syn_codons, self._ff_codons = self.__extract_syn_codons()

    def is_four_fold(self, codon):
        return codon in self._ff_codons

    def is_syn_codons(self, codon1: str, codon2: str):
        if not isinstance(codon1, str) or not isinstance(codon2, str):
            return False
        gc = self.codontable.forward_table
        return gc.get(codon1, "*") == gc.get(codon2, "*")
    
    def get_syn_number(self, cdn: str, pic: int):
        """return number of possible syn codons"""
        assert 0 <= pic <= 2, "pic must be 0-based and less than 3"
        return self._syn_codons.get((cdn, pic), 0)
        
    def get_mut_type(self, codon1: str, codon2: str, pic: int):
        """
        returned labels:
        - -1 - stopcodon loss or gain
        -  0 - usual sbs
        -  1 - synonimous sbs
        -  2 - fourfold sbs
        """
        assert codon1 != codon2, "codons must be different"
        assert 0 <= pic <= 2, "pic must be 0-based and less than 3"
        aa1 = self.codontable.forward_table.get(codon1, "*")
        aa2 = self.codontable.forward_table.get(codon2, "*")
        if aa1 == "*" or aa2 == "*":
            label = -1
        elif aa1 == aa2:
            label = 1
            if pic == 2 and self.is_four_fold(codon1):
                label = 2
        else:
            label = 0

        return label, aa1, aa2

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


def read_start_stop_codons(codontable: Union[NCBICodonTableDNA, int]):
    codontable = CodonAnnotation._prepare_codontable(codontable)
    return set(codontable.start_codons), set(codontable.stop_codons)


def calculate_mutspec(mutations: pd.DataFrame, freqs: Dict[str, float], label: str, use_context=True, use_proba=True):
    """
    mutations dataframe must contain 2 columns:
    - Mut
    - Label
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
