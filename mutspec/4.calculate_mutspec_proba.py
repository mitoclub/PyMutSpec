"""
pic is position in codon

TODO:
+ optimal states reading
- optimal mutations and spectra writing to files
- syn support
+ minimal sequence is gene (Part), not genome
+ 192 comp
- probability approach
- 
"""

import os
import sys
from collections import defaultdict
from datetime import datetime
from queue import Queue
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from utils import (
    iter_tree_edges, node_parent, profiler, get_farthest_leaf, 
    CodonAnnotation, GenomeStates,
    possible_codons, possible_sbs12, possible_sbs192
)
from custom_logging import logger


class MutSpec(CodonAnnotation, GenomeStates):
    nucl_order = ["A", "C", "G", "T"]
    # EPS = 1e-5
    
    def __init__(
            self, path_to_tree, path_to_anc, path_to_leaves, out_dir, 
            path_to_db="data/states.db", gcode=2, run=True, db_mode="db", 
            rewrite_db=False,
        ):
        CodonAnnotation.__init__(self, gcode)
        GenomeStates.__init__(self, path_to_anc, path_to_leaves, path_to_db, db_mode, rewrite_db)

        self.gcode = gcode
        self.MUT_LABELS = ["all", "ff"]  # TODO add syn
        self.tree = PhyloTree(path_to_tree, format=1)
        self.max_dist = get_farthest_leaf(self.tree)
        logger.info(
            f"tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
            f"max distance to leaf: {self.max_dist: .2f}"
        )
        if run:
            self.calc(path_to_anc, path_to_leaves, out_dir)

    def calc(self, path_to_anc, path_to_leaves, out_dir):
        mutations = self.extract_mutspec_from_tree(self.tree)
        # mutations, edge_mutspec12, edge_mutspec192, total_nucl_freqs = self.extract_mutspec_from_tree(self.tree)

        exit(0)
        os.makedirs(out_dir)
        print(f"Output directory {out_dir} created", file=sys.stderr)
        path_to_mutations = os.path.join(out_dir, "mutations.csv")
        path_to_nucl_freqs = os.path.join(out_dir, "nucl_freqs.csv")
        path_to_mutspec = os.path.join(out_dir, "mutspec{}_{}.csv")
        
        mutations.to_csv(path_to_mutations, index=None)
        total_nucl_freqs.to_csv(path_to_nucl_freqs, index=None)
        for label in self.MUT_LABELS:
            fp_mutspec12 = path_to_mutspec.format(12, label)
            fp_mutspec192 = path_to_mutspec.format(192, label)
            edge_mutspec12[label].to_csv(fp_mutspec12, index=None)
            edge_mutspec192[label].to_csv(fp_mutspec192, index=None)


    def extract_mutspec_from_tree(self, tree: PhyloTree):
        logger.info("Start mutation extraction from tree")
        edge_mutspec12 = defaultdict(list)  # all, syn, ff
        edge_mutspec192 = defaultdict(list)
        # full_tree_mutations = []  # TODO remove these lines, replace by online writing
        # total_nucl_freqs    = []  #
        for ref_node, alt_node in iter_tree_edges(tree):
            if ref_node.name not in self.nodes or alt_node.name not in self.nodes:
                continue
            logger.debug(f"extracting mutations from {ref_node.name} to {alt_node.name}")
            ref_genome = self.get_genome(ref_node.name)
            alt_genome  = self.get_genome(alt_node.name)
            _, dist_to_closest_leaf = ref_node.get_closest_leaf()
            evol_speed_coef = min(1, dist_to_closest_leaf / self.max_dist)

            genome_nucl_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
            genome_cxt_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
            genome_mutations = []
            genes_mutations = []
            for gene in ref_genome:
                ref_seq = ref_genome[gene]
                alt_seq = alt_genome[gene]

                gene_mut_df = self.extract_mutations(ref_seq, alt_seq)
                gene_mut_df["RefNode"] = ref_node.name
                gene_mut_df["AltNode"] = alt_node.name
                gene_mut_df["Gene"] = gene
                gene_mut_df["DistToLeaf"] = dist_to_closest_leaf
                gene_mut_df["EvolSpeedCoef"] = evol_speed_coef

                gene_nucl_freqs, gene_cxt_freqs = self.collect_state_freqs(ref_seq)
                for lbl in self.MUT_LABELS:
                    for nucl, freq in gene_nucl_freqs[lbl].items():
                        genome_nucl_freqs[lbl][nucl] += freq
                    for trinucl, freq in gene_cxt_freqs[lbl].items():
                        genome_cxt_freqs[lbl][trinucl] += freq

                if len(gene_mut_df) > 0:
                    genome_mutations.append(gene_mut_df)
                if len(gene_mut_df) > 200:
                    genes_mutations.append(gene_mut_df)  # TODO remove this line and 
            
            if len(genome_mutations) == 0:
                continue

            genome_mutations_df = pd.concat(genome_mutations)
            # full_tree_mutations.append(genome_mutations_df)
            
            # TODO add mutspec calculation

            cur_nucl_freqs = {"node": ref_node.name}
            for lbl in self.MUT_LABELS:
                for _nucl in "ACGT":
                    cur_nucl_freqs[f"{_nucl}_{lbl}"] = genome_nucl_freqs[lbl][_nucl]
                if lbl == "syn":
                    # TODO add syn support
                    raise NotImplementedError

                mutspec12 = self.calculate_mutspec12(genome_mutations_df, genome_nucl_freqs[lbl], label=lbl)
                mutspec12["RefNode"] = ref_node.name
                mutspec12["AltNode"] = alt_node.name
                mutspec192 = self.calculate_mutspec192(genome_mutations_df, genome_cxt_freqs[lbl], label=lbl)
                mutspec192["RefNode"] = ref_node.name
                mutspec192["AltNode"] = alt_node.name

                edge_mutspec12[lbl].append(mutspec12)
                edge_mutspec192[lbl].append(mutspec192)

            # total_nucl_freqs.append(cur_nucl_freqs)
        # if len(full_tree_mutations) == 10:
        #     break  # TODO remove line

        # full_tree_mutations_df = pd.concat(full_tree_mutations)
        # return full_tree_mutations_df
        # total_nucl_freqs_df = pd.DataFrame(total_nucl_freqs).drop_duplicates()  # TODO rewrite to normal optimal decision
        # edge_mutspec12_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec12.items()}
        # edge_mutspec192_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec192.items()}
        # return mutations_df, edge_mutspec12_df, edge_mutspec192_df, total_nucl_freqs_df


    def extract_mutations(self, g1: np.ndarray, g2: np.ndarray):
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
        - mut - dataframe of mutations
        - nucl_freqs - dict[lbl: dict[{ACGT}: int]] - nucleotide frequencies for all, syn and ff positions
        """
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        mutations = []
        for pos in range(1, n - 1):
            pic = pos % 3  # 0-based
            for cdn1, mut_cxt1, proba1 in self.sample_context_fast(pos, pic, g1):
                cdn1_str = "".join(cdn1)
                for cdn2, mut_cxt2, proba2 in self.sample_context_fast(pos, pic, g2):
                    cdn2_str = "".join(cdn2)

                    up_nuc1, up_nuc2 = mut_cxt1[0], mut_cxt2[0]
                    nuc1, nuc2 = mut_cxt1[1], mut_cxt2[1]
                    down_nuc1, down_nuc2 = mut_cxt1[2], mut_cxt2[2]

                    if nuc1 == nuc2 or up_nuc1 != up_nuc2 or down_nuc1 != down_nuc2:
                        continue
                    if sum([cdn1[_] == cdn2[_] for _ in range(3)]) != 2:
                        continue
                    
                    label, aa1, aa2 = self.get_mut_effect(cdn1_str, cdn2_str, pic)
                    sbs = {
                        "Mut": f"{up_nuc1}[{nuc1}>{nuc2}]{down_nuc1}",                        
                        "Effect": label,
                        "Pos": pos + 1,
                        "PosInCodon": pic + 1,
                        "RefCodon": cdn1_str,
                        "AltCodon": cdn2_str,
                        "RefAa": aa1,
                        "AltAa": aa2,
                        "ProbaRef": proba1,
                        "ProbaFull": proba1 * proba2,
                    }
                    mutations.append(sbs)

        mut_df = pd.DataFrame(mutations)
        return mut_df

    def collect_state_freqs(self, genome: np.ndarray):
        n = len(genome)
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        nucl_freqs = {lbl: defaultdict(np.float16) for lbl in ("all", "syn", "ff")}
        cxt_freqs = {lbl: defaultdict(np.float16) for lbl in ("all", "syn", "ff")}

        # nucl_freqs["all"] = {
        #     self.nucl_order[i]: fr for i, fr in enumerate(genome.sum(axis=1))
        # }
        for pos in range(1, n - 1):
            pic = pos % 3  # 0-based
            for cdn, cxt, proba in self.sample_context_fast(pos, pic, genome, 0.001):
                nuc = cxt[1]
                nucl_freqs["all"][nuc] += proba
                cxt_freqs["all"][cxt]  += proba

                _syn_scaler = self.get_syn_number(cdn, pic)
                if _syn_scaler > 0:
                    nucl_freqs["syn"][nuc] += proba * _syn_scaler
                    cxt_freqs["syn"][cxt]  += proba * _syn_scaler
                    if pic == 2 and self.is_four_fold(cdn):
                        nucl_freqs["ff"][nuc] += proba
                        cxt_freqs["ff"][cxt]  += proba

        return nucl_freqs, cxt_freqs

    @staticmethod
    def calculate_mutspec12(mut: pd.DataFrame, nucl_freqs, label: str):
        cols = ["Effect", "Mut"]
        for c in cols:
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

        mutspec = mut[mut["Effect"] >= label].Mut.value_counts().reset_index()
        mutspec.columns = ["Mut", "ObsFr"]

        mutspec_appendix = []
        unobserved_sbs = possible_sbs12.difference(mutspec.Mut.values)
        for usbs in unobserved_sbs:
            mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
        mutspec = pd.concat(
            [mutspec, pd.DataFrame(mutspec_appendix)],
            ignore_index=True
        )
        mutspec["RefNuc"] = mutspec.Mut.str.get(0)
        mutspec["Divisor"] = mutspec.RefNuc.map(nucl_freqs)
        mutspec["RawMutSpec"] = mutspec.ObsFr / mutspec.Divisor
        mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
        mutspec.drop("RefNuc", axis=1, inplace=True)
        return mutspec

    @staticmethod
    def calculate_mutspec192(mut: pd.DataFrame, codon_freqs, label: str):
        cols = ["Effect", "MutExt"]
        for c in cols:
            assert c in mut.columns, f"Column {c} is not in mut df"

        available_labels = {"syn", "ff", "all"}
        if isinstance(label, str):
            label = label.lower()
            if label not in available_labels:
                raise ValueError(f"pass the appropriate label: {available_labels}")
            if label == "syn":
                label = 1
            elif label == "ff":
                label = 2
            elif label == "all":
                label = 0
        else:
            raise ValueError(f"pass the appropriate label: {available_labels}")

        mutspec = mut[mut["Effect"] >= label].MutExt.value_counts().reset_index()
        mutspec.columns = ["Mut", "ObsFr"]

        mutspec_appendix = []
        unobserved_sbs = possible_sbs192.difference(mutspec.Mut.values)
        for usbs in unobserved_sbs:
            mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
        mutspec = pd.concat(
            [mutspec, pd.DataFrame(mutspec_appendix)],
            ignore_index=True
        )
        mutspec["Context"] = mutspec.Mut.str.get(0) + mutspec.Mut.str.get(2) + mutspec.Mut.str.get(-1)
        mutspec["Divisor"] = mutspec.Context.map(codon_freqs)
        mutspec["RawMutSpec"] = (mutspec.ObsFr / mutspec.Divisor).fillna(0)
        mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
        mutspec.drop("Context", axis=1, inplace=True)
        return mutspec

    def sample_context_fast(self, pos, pic, genome: np.ndarray, cutoff=0.01):
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

            codon = tuple(self.nucl_order[_] for _ in (i, j, k))
            if pic == 0:
                mut_context = tuple(self.nucl_order[_] for _ in (m, i, j))
                full_proba = probas[m, k, j, i]
            elif pic == 2:
                mut_context = tuple(self.nucl_order[_] for _ in (j, k, m))
                full_proba = probas[m, k, j, i]
            elif pic == 1:
                mut_context = codon
                full_proba = probas[k, j, i]
            
            yield codon, mut_context, full_proba


def main():
    path_to_tree =   "./data/example_birds/anc_kg.treefile"
    path_to_states = "./data/example_birds/anc_kg_states_birds.tsv"
    path_to_leaves = "./data/example_birds/leaves_birds_states.tsv"
    out_dir = "./data/processed/birds"
    out_dir = out_dir + "_" + datetime.now().strftime("%d-%m-%y-%H-%M-%S")
    MutSpec(path_to_tree, path_to_states, path_to_leaves, out_dir, run=True)


if __name__ == "__main__":
    main()
