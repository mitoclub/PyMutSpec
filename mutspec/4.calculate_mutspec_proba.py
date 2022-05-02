"""
TODO:
- optimal states reading
- optimal mutations and spectra writing to files
- syn support
+ minimal sequence is gene (Part), not genome
+ 192 comp
- probability approach
- 
"""

import logging
import logging.config
import os
import sys
from collections import defaultdict
from datetime import datetime
from queue import Queue
from typing import Iterable

import yaml
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from utils import (extract_syn_codons, is_syn_codons, node_parent, profiler,
                   possible_codons, possible_sbs12, possible_sbs192)

EPS = 1e-5
DEFAULT_PATH_TO_LOGCONF = "./configs/log_settings.yaml"


with open(os.environ.get("LOG_CONFIG", DEFAULT_PATH_TO_LOGCONF), "rb") as file:
    log_config = yaml.safe_load(file)
    logging.config.dictConfig(log_config)
    del log_config

logger = logging.getLogger('MutSpecCalc')


class MutSpec:
    def __init__(
            self, 
            path_to_tree,
            path_to_states,
            path_to_leaves,
            out_dir,
            gcode=2,
        ):
        self.MUT_LABELS = ["all", "ff"]  # TODO add syn
        self.nucl_order = ["A", "C", "G", "T"]
        self.gcode = gcode
        self.codontable = CodonTable.unambiguous_dna_by_id[gcode]
        self.ff_codons, self.syn_codons = extract_syn_codons(self.codontable)

        tree, anc, leaves = self.read_data(path_to_tree, path_to_states, path_to_leaves)
        self.node2genome = self.precalc_node2genome(anc, leaves)
        self.nodes = set(self.node2genome.keys())
        mutations = self.extract_mutspec_from_tree(tree)
        # mutations, edge_mutspec12, edge_mutspec192, total_nucl_freqs = self.extract_mutspec_from_tree(tree)

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

    @staticmethod
    def precalc_node2genome(anc: pd.DataFrame, leaves: pd.DataFrame) -> dict:
        node2genome = defaultdict(dict)
        for states in [anc, leaves]:
            logger.debug("Sorting states")
            states = states.sort_values(["Node", "Part", "Site"])
            logger.debug("States sorted")
            gr = states.groupby(["Node", "Part"])
            for (node, part), gene_pos_ids in gr.groups.items():
                gene_df = states.loc[gene_pos_ids]
                gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
                node2genome[node][part] = gene_states
        logger.info("node2genome mapping builded")
        return node2genome
    
    def get_genome(self, node: str):
        return self.node2genome[node]

    @staticmethod
    def read_data(path_to_tree, path_to_states, path_to_leaves, states_dtype=np.float16):
        tree = PhyloTree(path_to_tree, format=1)
        logger.info(f"tree loaded, number of leaf nodes: {len(tree)}, total number of nodes: {len(tree.get_cached_content())}")
        dtypes = {
            "p_A": states_dtype, 
            "p_C": states_dtype, 
            "p_G": states_dtype, 
            "p_T": states_dtype,
            "Site": np.int32,
            "Part": np.int8,
        }
        usecols = ["Node", "Part", "Site", "p_A", "p_C", "p_G", "p_T"]
        anc = pd.read_csv(path_to_states, sep="\t", comment='#', usecols=usecols, dtype=dtypes)
        logger.info(f"anc states loaded, shape = {anc.shape}")
        leaves = pd.read_csv(path_to_leaves, sep="\t", usecols=usecols, dtype=dtypes)
        logger.info(f"leaves states loaded, shape = {leaves.shape}")
        aln_sizes = leaves.groupby("Node").apply(len)
        assert aln_sizes.nunique() == 1, "uncomplete leaves state table"
        return tree, anc, leaves

    def extract_mutspec_from_tree(self, tree):
        logger.info("Start extracting mutspec from tree")
        discovered_nodes = set()
        discovered_nodes.add(tree.name)
        Q = Queue()
        Q.put(tree)

        edge_mutspec12 = defaultdict(list)  # all, syn, ff
        edge_mutspec192 = defaultdict(list)
        full_tree_mutations = []
        total_nucl_freqs = []
        while not Q.empty():
            cur_node = Q.get()
            for child in cur_node.children:
                Q.put(child)

            if cur_node.name not in discovered_nodes:
                discovered_nodes.add(cur_node.name)
                if cur_node.name not in self.nodes:
                    continue
                parent_node = node_parent(cur_node)

                # main process starts here
                logger.debug(f"extracting mutations from {parent_node.name} to {cur_node.name}")
                parent_gene = self.get_genome(parent_node.name)
                child_gene  = self.get_genome(cur_node.name)
                _, dist_to_closest_leaf = parent_node.get_closest_leaf()

                genome_nucl_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
                genome_cxt_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
                genome_mutations = []
                genes_mutations = []
                for gene in parent_gene:
                    gene_ref = parent_gene[gene]
                    gene_alt = child_gene[gene]
                    gene_mut_df = self.extract_mutations(
                        gene_ref, gene_alt,
                        parent_node.name, cur_node.name, gene,
                    )
                    gene_mut_df["DistToLeaf"] = dist_to_closest_leaf
                    gene_nucl_freqs, gene_cxt_freqs = self.collect_freqs(gene_ref)
                    for lbl in self.MUT_LABELS:
                        for nucl, freq in gene_nucl_freqs[lbl].items():
                            genome_nucl_freqs[lbl][nucl] += freq
                        for trinucl, freq in gene_cxt_freqs[lbl].items():
                            genome_cxt_freqs[lbl][trinucl] += freq

                    if len(gene_mut_df) > 0:
                        genome_mutations.append(gene_mut_df)
                    if len(gene_mut_df) > 200:
                        genes_mutations.append(gene_mut_df)
                
                if len(genome_mutations) == 0:
                    continue

                genome_mutations_df = pd.concat(genome_mutations)
                full_tree_mutations.append(genome_mutations_df)
                
                # TODO add mutspec calculation

                cur_nucl_freqs = {"node": parent_node.name}
                for lbl in self.MUT_LABELS:
                    for _nucl in "ACGT":
                        cur_nucl_freqs[f"{_nucl}_{lbl}"] = genome_nucl_freqs[lbl][_nucl]
                    if lbl == "syn":
                        raise NotImplementedError

                    mutspec12 = self.calculate_mutspec12(genome_mutations_df, genome_nucl_freqs[lbl], label=lbl)
                    mutspec12["RefNode"] = parent_node.name
                    mutspec12["AltNode"] = cur_node.name
                    mutspec192 = self.calculate_mutspec192(genome_mutations_df, genome_cxt_freqs[lbl], label=lbl)
                    mutspec192["RefNode"] = parent_node.name
                    mutspec192["AltNode"] = cur_node.name

                    edge_mutspec12[lbl].append(mutspec12)
                    edge_mutspec192[lbl].append(mutspec192)

                total_nucl_freqs.append(cur_nucl_freqs)
            # if len(full_tree_mutations) == 10:
            #     break  # TODO remove line

        full_tree_mutations_df = pd.concat(full_tree_mutations)
        return full_tree_mutations_df
        # total_nucl_freqs_df = pd.DataFrame(total_nucl_freqs).drop_duplicates()  # TODO rewrite to normal optimal decision
        # edge_mutspec12_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec12.items()}
        # edge_mutspec192_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec192.items()}
        # return mutations_df, edge_mutspec12_df, edge_mutspec192_df, total_nucl_freqs_df


    def extract_mutations(
            self, 
            g1: np.ndarray, g2: np.ndarray, 
            name1: str, name2: str, 
            gene: str,
            collect_freqs=True, context=False,
        ):
        """
        Extract alterations of g2 comparing to g1

        params:
        - g1 - reference sequence (parent node)
        - g2 - alternative sequence (child node)
        - name1 - node name of ref
        - name2 - node name of alt
        - collect_freqs - 
        - context - TODO

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
            pos_in_codon = pos % 3  # 0-based
            for cdn1, mut_cxt1, proba1 in self.sample_context_fast(pos, pos_in_codon, g1):
                cdn1_str = "".join(cdn1)
                for cdn2, mut_cxt2, proba2 in self.sample_context_fast(pos, pos_in_codon, g2):
                    cdn2_str = "".join(cdn2)

                    nuc1, nuc2 = mut_cxt1[1], mut_cxt2[1]
                    up_nuc1, up_nuc2 = mut_cxt1[0], mut_cxt2[0]
                    down_nuc1, down_nuc2 = mut_cxt1[2], mut_cxt2[2]

                    if nuc1 == nuc2:
                        continue
                    if up_nuc1 != up_nuc2 or down_nuc1 != down_nuc2:
                        continue
                    if sum([cdn1[_] == cdn2[_] for _ in range(3)]) != 2:
                        continue
                    
                    label, aa1, aa2 = self.get_mut_label(cdn1_str, cdn2_str, pos_in_codon)
                    sbs = {
                        "RefNode": name1,
                        "AltNode": name2,
                        "Gene": gene,
                        "Mut": f"{up_nuc1}[{nuc1}>{nuc2}]{down_nuc1}",                        
                        "Effect": label,
                        "Pos": pos + 1,
                        "PosInCodon": pos_in_codon + 1,
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

    def collect_freqs(self, genome: np.ndarray):
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
                if pic == 2 and cdn in self.ff_codons:
                    nucl_freqs["ff"][nuc] += proba
                    cxt_freqs["ff"][cxt]  += proba
                if (cdn, pic) in self.syn_codons:
                    nucl_freqs["syn"][nuc] += proba * self.syn_codons[(cdn, pic)]
                    cxt_freqs["syn"][cxt]  += proba * self.syn_codons[(cdn, pic)]

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

    def sample_context(self, pos, pos_in_codon, genome: np.ndarray, cutoff=0.01):
        nuc_cutoff = cutoff * 5
        codon_states = genome[pos - pos_in_codon: pos - pos_in_codon + 3]        
        extra_codon_states = genome[pos + pos_in_codon - 1]  # doesn't mean if pos_in_codon == 1
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
                    
                    if pos_in_codon != 1:
                        for m, p4 in enumerate(extra_codon_states):
                            if p4 < nuc_cutoff:
                                continue
                            full_proba = codon_proba * p4
                            if full_proba < cutoff:
                                continue
                            if pos_in_codon == 0:
                                mut_context = tuple(self.nucl_order[_] for _ in (m, i, j))
                            elif pos_in_codon == 2:
                                mut_context = tuple(self.nucl_order[_] for _ in (j, k, m))
                            yield codon, mut_context, full_proba
                    else:
                        yield codon, codon, codon_proba

    def sample_context_fast(self, pos, pos_in_codon, genome: np.ndarray, cutoff=0.01):
        codon_states = genome[pos - pos_in_codon: pos - pos_in_codon + 3]        
        extra_codon_states = genome[pos + pos_in_codon - 1]  # doesn't mean if pos_in_codon == 1
        # gaps are not appropriate
        if np.any(codon_states.sum(1) == 0) or extra_codon_states.sum() == 0:
            return

        a, b, c = codon_states
        probas = a * b[:, None] * c[:, None, None]
        _ii = 0  # increment index if there are 4th context nucl
        if pos_in_codon != 1:
            probas = probas * extra_codon_states[:, None, None, None]
            _ii = 1

        indexes = np.where(probas > cutoff)
        for idx in range(len(indexes[0])):
            i, j, k = indexes[2+_ii][idx], indexes[1+_ii][idx], indexes[0+_ii][idx]
            m = indexes[0][idx]

            codon = tuple(self.nucl_order[_] for _ in (i, j, k))
            if pos_in_codon == 0:
                mut_context = tuple(self.nucl_order[_] for _ in (m, i, j))
                full_proba = probas[m, k, j, i]
            elif pos_in_codon == 2:
                mut_context = tuple(self.nucl_order[_] for _ in (j, k, m))
                full_proba = probas[m, k, j, i]
            elif pos_in_codon == 1:
                mut_context = codon
                full_proba = probas[k, j, i]
            
            yield codon, mut_context, full_proba

    def is_four_fold(self, codon):
        return codon in self.ff_codons

    def is_syn(self, codon1, codon2):
        return self.codontable.forward_table[codon1] == self.codontable.forward_table[codon2]
        
    def get_mut_label(self, codon1: str, codon2: str, pos_in_codon: int):
        """
        returned labels:
        - -1 - error sbs (stopcodon alterations)
        -  0 - usual sbs
        -  1 - synonimous sbs
        -  2 - fourfold sbs
        """
        assert codon1 != codon2, "codons must be different"
        aa1 = self.codontable.forward_table.get(codon1, "*")
        aa2 = self.codontable.forward_table.get(codon2, "*")
        if aa1 == "*" or aa2 == "*":
            label = -1
        elif aa1 == aa2:
            label = 1
            if pos_in_codon == 2 and self.is_four_fold(codon1):
                label = 2
        else:
            label = 0

        return label, aa1, aa2

def main():
    path_to_tree =   "./data/interim/iqtree_runs/brun3/anc_kg.treefile"
    path_to_states = "./data/interim/anc_kg_states_birds.tsv"
    path_to_leaves = "./data/interim/leaves_birds_states.tsv"
    out_dir = "./data/processed/birds"
    out_dir = out_dir + "_" + datetime.now().strftime("%d-%m-%y-%H-%M-%S")
    MutSpec(path_to_tree, path_to_states, path_to_leaves, out_dir)


if __name__ == "__main__":
    main()
