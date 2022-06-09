"""
pic is position in codon
"""

import os
import sys
from collections import defaultdict
from datetime import datetime
from queue import Queue
from typing import Dict, Iterable

import click
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from utils import (
    iter_tree_edges, profiler, calculate_mutspec,
    CodonAnnotation, GenomeStates, possible_codons
)
from custom_logging import logger


class MutSpec(CodonAnnotation, GenomeStates):    
    def __init__(
            self, path_to_tree, path_to_states, out_dir, 
            path_to_db="data/states.db", gcode=2, run=True, db_mode="dict", 
            rewrite_db=None, proba_cutoff=0.01, proba_mode=False,
        ):
        for path in path_to_states + [path_to_tree]:
            if not os.path.exists(path):
                raise ValueError(f"Path doesn't exist: {path}")
        if os.path.exists(out_dir):
            raise ValueError(f"Out directory path exist: {path}")

        CodonAnnotation.__init__(self, gcode)
        GenomeStates.__init__(
            self, path_to_states, path_to_db, db_mode, rewrite_db, proba_mode
        )
        self.gcode = gcode
        logger.info(f"Using gencode {gcode}")
        self.proba_cutoff = proba_cutoff
        self.MUT_LABELS = ["all", "syn", "ff"]
        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        # self.max_dist = self.fp_format(get_farthest_leaf(self.tree))
        logger.info(
            f"tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
            # f"max distance to leaf: {self.max_dist: .2f}"
        )
        if run:
            self.open_handles(out_dir)
            self.extract_mutspec_from_tree()
            self.close_handles()

    def open_handles(self, out_dir):
        os.makedirs(out_dir)
        logger.info(f"Output directory '{out_dir}' created")
        path_to_mutations  = os.path.join(out_dir, "mutations.tsv")
        path_to_nucl_freqs = os.path.join(out_dir, "freqs.tsv")
        path_to_mutspec12  = os.path.join(out_dir, "mutspec12.tsv")
        path_to_mutspec192 = os.path.join(out_dir, "mutspec192.tsv")
        path_to_mutspec_genes12  = os.path.join(out_dir, "mutspec12genes.tsv")
        path_to_mutspec_genes192 = os.path.join(out_dir, "mutspec192genes.tsv")
        self.handle = dict()
        self.handle["mut"]  = open(path_to_mutations, "w")
        self.handle["freq"] = open(path_to_nucl_freqs, "w")
        self.handle["ms12"] = open(path_to_mutspec12, "w")
        self.handle["ms192"] = open(path_to_mutspec192, "w")
        self.handle["ms12s"] = open(path_to_mutspec_genes12, "w")
        self.handle["ms192g"] = open(path_to_mutspec_genes192, "w")
        logger.info("Handles opened")

    def close_handles(self):
        for x in self.handle:
            self.handle[x].close()
        logger.info("Handles closed")

    # @profiler
    def extract_mutspec_from_tree(self):
        logger.info("Start mutation extraction from tree")
        add_header = defaultdict(lambda: True)
        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree), 1):
            if ref_node.name not in self.nodes or alt_node.name not in self.nodes:
                continue
            logger.debug(f"Extracting mutations from {ref_node.name} to {alt_node.name}")

            # get genomes from storage
            ref_genome  = self.get_genome(ref_node.name)
            alt_genome  = self.get_genome(alt_node.name)

            # calculate phylogenetic uncertainty correction
            # _, dist_to_closest_leaf = ref_node.get_closest_leaf()
            # dist_to_closest_leaf = self.fp_format(dist_to_closest_leaf)
            # evol_speed_coef = self.fp_format(1 - min(1, dist_to_closest_leaf / self.max_dist))

            genome_nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_cxt_freqs  = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_mutations = []
            for gene in ref_genome:
                ref_seq = ref_genome[gene]
                alt_seq = alt_genome[gene]
                gene = np.int16(gene)
                
                # extract mutations and put in order columns
                gene_mut_df = self.extract_mutations_simple(ref_seq, alt_seq)
                # gene_mut_df["DistToClosestLeaf"] = dist_to_closest_leaf
                # gene_mut_df["EvolSpeedCoef"] = evol_speed_coef
                # gene_mut_df["ProbaFull"] = gene_mut_df["EvolSpeedCoef"] * gene_mut_df["ProbaMut"]
                gene_mut_df["RefNode"] = ref_node.name
                gene_mut_df["AltNode"] = alt_node.name
                gene_mut_df["Gene"] = gene

                # collect state frequencies
                gene_nucl_freqs, gene_cxt_freqs = self.collect_obs_mut_freqs(ref_seq)

                # dump state frequencies
                self.dump_freqs(
                    gene_nucl_freqs, gene_cxt_freqs, ref_node.name, 
                    gene, self.handle["freq"], add_header["freqs"],
                )
                add_header["freqs"] = False

                # summarize state frequencies over genome
                for lbl in self.MUT_LABELS:
                    for nucl, freq in gene_nucl_freqs[lbl].items():
                        genome_nucl_freqs[lbl][nucl] += freq
                    for trinucl, freq in gene_cxt_freqs[lbl].items():
                        genome_cxt_freqs[lbl][trinucl] += freq
                
                # calculate gene mutational spectra for all labels
                if len(gene_mut_df) > 0:
                    genome_mutations.append(gene_mut_df)
                    for lbl in self.MUT_LABELS:
                        mutspec12 = calculate_mutspec(gene_mut_df, gene_nucl_freqs[lbl], label=lbl, use_context=False, use_proba=False)
                        mutspec12["RefNode"] = ref_node.name
                        mutspec12["AltNode"] = alt_node.name
                        mutspec12["Label"] = lbl
                        mutspec12["Gene"]  = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec12, self.handle["ms12s"], add_header["ms12g"])
                        add_header["ms12g"] = False

                if len(gene_mut_df) > 100:
                    for lbl in self.MUT_LABELS:
                        mutspec192 = calculate_mutspec(gene_mut_df, gene_cxt_freqs[lbl], label=lbl, use_context=True, use_proba=False)
                        mutspec192["RefNode"] = ref_node.name
                        mutspec192["AltNode"] = alt_node.name
                        mutspec192["Label"] = lbl
                        mutspec192["Gene"] = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec192, self.handle["ms192g"], add_header["ms192g"])
                        add_header["ms192g"] = False

            if len(genome_mutations) == 0:
                continue

            genome_mutations_df = pd.concat(genome_mutations)
            del genome_mutations

            # dump mutations
            self.dump_table(genome_mutations_df, self.handle["mut"], add_header["mut"])
            add_header["mut"] = False
            
            # calculate full genome mutational spectra for all labels
            for lbl in self.MUT_LABELS:
                mutspec12 = calculate_mutspec(genome_mutations_df, genome_nucl_freqs[lbl], label=lbl, use_context=False, use_proba=False)
                mutspec12["RefNode"] = ref_node.name
                mutspec12["AltNode"] = alt_node.name
                mutspec12["Label"] = lbl
                mutspec192 = calculate_mutspec(genome_mutations_df, genome_cxt_freqs[lbl], label=lbl, use_context=True, use_proba=False)
                mutspec192["RefNode"] = ref_node.name
                mutspec192["AltNode"] = alt_node.name
                mutspec192["Label"] = lbl

                # Dump genome mutspecs
                self.dump_table(mutspec12,  self.handle["ms12"],  add_header["ms"])
                self.dump_table(mutspec192, self.handle["ms192"], add_header["ms"])
                add_header["ms"] = False

            # if ei == 3:
            #     break  # TODO remove lines

            if ei % 100 == 0:
                logger.info(f"Processed {ei} tree edges")

        logger.info(f"Processed {ei} tree edges")
        logger.info("MutSpec extraction done")

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
    

# @click.option("--out", "out_dir", required=True, type=click.Path(writable=True), help="")
# @click.command("MutSpec calculator", help="")
# @click.option("--tree", "path_to_tree", required=True, type=click.Path(True), help="")
# @click.option("--anc", "path_to_anc", required=True, type=click.Path(True), help="")
# @click.option("--leaves", "path_to_leaves", required=True, type=click.Path(True), help="path to leaves table of custom format")
def main(path_to_tree, path_to_anc, path_to_leaves, out_dir_lbl=""):
    out_dir = "./data/processed/nematoda"
    out_dir_lbl = datetime.now().strftime("%d-%m-%y-%H-%M-%S") + "_" + out_dir_lbl
    out_dir = os.path.join(out_dir, out_dir_lbl)
    MutSpec(path_to_tree, [path_to_anc, path_to_leaves], out_dir, run=True, gcode=5)


if __name__ == "__main__":
    path_to_tree =   "./data/example_nematoda/anc.treefile.rooted"
    path_to_leaves = "./data/example_nematoda/leaves_states_nematoda.tsv"
    path_to_anc = "./data/example_nematoda/nematoda_anc_mf/genes_states.tsv"
    main(path_to_tree, path_to_anc, path_to_leaves, "anc_mf")
