"""
pic is position in codon

TODO:
+ optimal states reading
+ optimal mutations and spectra writing to files
+ syn support
+ minimal sequence is gene (Part), not genome
+ 192 comp
+ probability approach
- 
"""

import os
import sys
from collections import defaultdict
from datetime import datetime
from queue import Queue
from typing import Dict, Iterable

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from mutspec.utils import (
    iter_tree_edges, profiler, get_farthest_leaf, calculate_mutspec,
    CodonAnnotation, GenomeStates, possible_codons
)
from _custom_logging import logger


class MutSpec(CodonAnnotation, GenomeStates):
    nucl_order = ["A", "C", "G", "T"]
    # EPS = 1e-5
    
    def __init__(
            self, path_to_tree, path_to_states, out_dir, 
            path_to_db="data/states.db", gcode=2, run=True, db_mode="dict", 
            rewrite_db=None, proba_cutoff=0.01,
        ):
        for path in path_to_states + [path_to_tree]:
            if not os.path.exists(path):
                raise ValueError(f"Path doesn't exist: '{path}'")
        # if os.path.exists(out_dir):
        #     raise ValueError(f"Out directory path exist: '{path}'")

        CodonAnnotation.__init__(self, gcode)
        GenomeStates.__init__(self, path_to_states, path_to_db, db_mode, rewrite_db, True)

        self.gcode = gcode
        logger.info(f"Using gencode {gcode}")
        self.proba_cutoff = proba_cutoff
        self.MUT_LABELS = ["all", "syn", "ff"]
        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        self.max_dist = self.fp_format(get_farthest_leaf(self.tree))
        logger.info(
            f"tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
            f"max distance to leaf: {self.max_dist: .2f}"
        )
        if run:
            self.calc(out_dir)

    def calc(self, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        logger.info(f"Output directory '{out_dir}' created")
        path_to_mutations = os.path.join(out_dir, "mutations.csv")
        path_to_nucl_freqs = os.path.join(out_dir, "freqs.csv")
        path_to_mutspec12 = os.path.join(out_dir, "mutspec12.csv")
        path_to_mutspec192 = os.path.join(out_dir, "mutspec192.csv")
        path_to_mutspec_genes12 = os.path.join(out_dir, "mutspec_genes12.csv")
        path_to_mutspec_genes192 = os.path.join(out_dir, "mutspec_genes192.csv")

        self.handle_mut = open(path_to_mutations, "w")
        self.handle_freq = open(path_to_nucl_freqs, "w")
        self.handle_mutspec12 = open(path_to_mutspec12, "w")
        self.handle_mutspec192 = open(path_to_mutspec192, "w")
        self.handle_mutspec_genes12 = open(path_to_mutspec_genes12, "w")
        self.handle_mutspec_genes192 = open(path_to_mutspec_genes192, "w")

        self.extract_mutspec_from_tree()

        self.handle_mut.close()
        self.handle_freq.close()
        self.handle_mutspec12.close()
        self.handle_mutspec192.close()
        self.handle_mutspec_genes12.close()
        self.handle_mutspec_genes192.close()

    # @profiler
    def extract_mutspec_from_tree(self):
        logger.info("Start mutation extraction from tree")
        add_header = defaultdict(lambda: True)
        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree)):
            # if alt_node.name.startswith("Node"):
            #     continue  # TODO drop line
            if ref_node.name not in self.nodes or alt_node.name not in self.nodes:
                logger.debug(f"Pass edge '{ref_node.name}'-'{alt_node.name}' due to absence of genome")
                continue
            logger.debug(f"extracting mutations from {ref_node.name} to {alt_node.name}")

            # get genomes from storage
            ref_genome = self.get_genome(ref_node.name)
            alt_genome  = self.get_genome(alt_node.name)

            # calculate phylogenetic uncertainty correction
            _, dist_to_closest_leaf = ref_node.get_closest_leaf()
            dist_to_closest_leaf = self.fp_format(dist_to_closest_leaf)
            evol_speed_coef = self.fp_format(1 - min(1, dist_to_closest_leaf / self.max_dist))

            genome_nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_cxt_freqs = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_mutations = []
            for gene in ref_genome:
                ref_seq = ref_genome[gene]
                alt_seq = alt_genome[gene]
                # gene = np.int16(gene)
                
                # extract mutations and put in order columns
                gene_mut_df = self.extract_mutations(ref_seq, alt_seq)
                if gene_mut_df.shape[0] == 0:
                    continue
                # gene_mut_df["DistToClosestLeaf"] = dist_to_closest_leaf
                # gene_mut_df["EvolSpeedCoef"] = evol_speed_coef
                # gene_mut_df["ProbaFull"] = gene_mut_df["EvolSpeedCoef"] * gene_mut_df["ProbaMut"]
                gene_mut_df["ProbaFull"] = gene_mut_df["ProbaMut"]
                gene_mut_df["RefNode"] = ref_node.name
                gene_mut_df["AltNode"] = alt_node.name
                gene_mut_df["Gene"] = gene

                # collect state frequencies
                gene_nucl_freqs, gene_cxt_freqs = self.collect_state_freqs(ref_seq, evol_speed_coef)

                # dump state frequencies
                self.dump_freqs(
                    gene_nucl_freqs, gene_cxt_freqs, ref_node.name, 
                    gene, self.handle_freq, add_header["freqs"],
                )
                add_header["freqs"] = False
                # summarize state frequencies over genome
                for lbl in self.MUT_LABELS:
                    for nucl, freq in gene_nucl_freqs[lbl].items():
                        genome_nucl_freqs[lbl][nucl] += freq
                    for trinucl, freq in gene_cxt_freqs[lbl].items():
                        genome_cxt_freqs[lbl][trinucl] += freq

                if len(gene_mut_df) > 0:
                    genome_mutations.append(gene_mut_df)
                
                # calculate gene mutational spectra for all labels
                if len(gene_mut_df) > 0:
                    for lbl in self.MUT_LABELS:
                        mutspec12 = calculate_mutspec(gene_mut_df, gene_nucl_freqs[lbl], label=lbl, use_context=False)
                        mutspec12["RefNode"] = ref_node.name
                        mutspec12["AltNode"] = alt_node.name
                        mutspec12["Label"] = lbl
                        mutspec12["Gene"] = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec12,  self.handle_mutspec_genes12,  add_header["ms12g"])
                        add_header["ms12g"] = False

                if len(gene_mut_df) > 100:  # TODO
                    for lbl in self.MUT_LABELS:
                        mutspec192 = calculate_mutspec(gene_mut_df, gene_cxt_freqs[lbl], label=lbl, use_context=True)
                        mutspec192["RefNode"] = ref_node.name
                        mutspec192["AltNode"] = alt_node.name
                        mutspec192["Label"] = lbl
                        mutspec192["Gene"] = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec192, self.handle_mutspec_genes192, add_header["ms192g"])
                        add_header["ms192g"] = False
            
            if len(genome_mutations) == 0:
                continue

            genome_mutations_df = pd.concat(genome_mutations)
            del genome_mutations

            # dump mutations
            self.dump_table(genome_mutations_df, self.handle_mut, add_header["mut"])
            add_header["mut"] = False
            
            # calculate full genome mutational spectra for all labels
            for lbl in self.MUT_LABELS:
                mutspec12 = calculate_mutspec(genome_mutations_df, genome_nucl_freqs[lbl], label=lbl, use_context=False)
                mutspec12["RefNode"] = ref_node.name
                mutspec12["AltNode"] = alt_node.name
                mutspec12["Label"] = lbl
                mutspec192 = calculate_mutspec(genome_mutations_df, genome_cxt_freqs[lbl], label=lbl, use_context=True)
                mutspec192["RefNode"] = ref_node.name
                mutspec192["AltNode"] = alt_node.name
                mutspec192["Label"] = lbl

                # Dump genome mutspecs
                self.dump_table(mutspec12,  self.handle_mutspec12,  add_header["ms"])
                self.dump_table(mutspec192, self.handle_mutspec192, add_header["ms"])
                add_header["ms"] = False

            if ei % 100 == 0:
                logger.info(f"Processed {ei} tree edges")
            # if ei > 10:
            #     break  # TODO remove lines

    def extract_mutations(self, g1: np.ndarray, g2: np.ndarray):
        """
        Extract alterations of g2 comparing to g1

        Arguments
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
            for cdn1, mut_cxt1, proba1 in self.sample_context_fast(pos, pic, g1, self.proba_cutoff):
                cdn1_str = "".join(cdn1)
                for cdn2, mut_cxt2, proba2 in self.sample_context_fast(pos, pic, g2, self.proba_cutoff):
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
                    }
                    mutations.append(sbs)

        mut_df = pd.DataFrame(mutations)
        return mut_df

    def collect_state_freqs(self, genome: np.ndarray, evol_speed_coef: float,  proba_cutoff=0.001):
        n = len(genome)
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        assert 0 <= evol_speed_coef <= 1, "Evol coefficient must be between 0 and 1"

        nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in ("all", "syn", "ff")}
        cxt_freqs = {lbl: defaultdict(self.fp_format) for lbl in ("all", "syn", "ff")}

        for pos in range(1, n - 1):
            pic = pos % 3  # 0-based
            for cdn, cxt, proba in self.sample_context_fast(pos, pic, genome, proba_cutoff):
                cdn_str = "".join(cdn)
                proba *= evol_speed_coef  # adjusted proba
                nuc = cxt[1]
                cxt_str = "".join(cxt)

                nucl_freqs["all"][nuc] += proba
                cxt_freqs["all"][cxt_str]  += proba

                _syn_scaler = self.get_syn_number(cdn_str, pic)
                if _syn_scaler > 0:
                    nucl_freqs["syn"][nuc] += proba * _syn_scaler
                    cxt_freqs["syn"][cxt_str]  += proba * _syn_scaler
                    if pic == 2 and self.is_four_fold(cdn_str):
                        nucl_freqs["ff"][nuc] += proba
                        cxt_freqs["ff"][cxt_str]  += proba

        return nucl_freqs, cxt_freqs

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

    @staticmethod
    def dump_table(df: pd.DataFrame, handle, header=False):
        if header:
            handle.write("\t".join(df.columns) + "\n")
        handle.write(df.to_csv(sep="\t", index=None, header=None))
    
    def dump_freqs(self, gene_nuc_freqs, gene_cxt_freqs, node, gene, handle, header=False):
        if header:
            nucls  = "\t".join(self.nucl_order)
            codons = "\t".join(possible_codons)
            header = f"Node\tGene\tLabel\t{nucls}\t{codons}\n"
            handle.write(header)
        
        for lbl in self.MUT_LABELS:
            nucls  = "\t".join([str(gene_nuc_freqs[lbl][x]) for x in self.nucl_order])
            codons = "\t".join([str(gene_cxt_freqs[lbl][x]) for x in possible_codons])
            row = f"{node}\t{gene}\t{lbl}\t{nucls}\t{codons}\n"
            handle.write(row)
    

def main():
    path_to_tree =   "./data/example_nematoda/anc.treefile.rooted"
    path_to_anc = "./data/example_nematoda/genes_states.pastml.tsv"
    # path_to_leaves = "./data/example_nematoda/leaves_states_nematoda.tsv"
    out_dir = "./data/processed/nematoda/pastml"
    # out_dir = os.path.join(out_dir, datetime.now().strftime("%d-%m-%y-%H-%M-%S") + "_proba")
    MutSpec(
        path_to_tree, [path_to_anc], 
        out_dir, run=True, db_mode="dict", gcode=5, 
    )


if __name__ == "__main__":
    main()
