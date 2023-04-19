#!/usr/bin/env python3
"""
pic is position in codon

"""

import os
import sys
from collections import defaultdict
from shutil import rmtree
from typing import Iterable, Union

import click
import numpy as np
import pandas as pd
from ete3 import PhyloTree

from pymutspec.annotation import (CodonAnnotation, calculate_mutspec,
                                      get_farthest_leaf, iter_tree_edges,
                                      lbl2lbl_id)
from pymutspec.constants import possible_sbs12, possible_sbs192
from pymutspec.io import GenesStates
from pymutspec.utils import load_logger

logger = None


class MutSpec(CodonAnnotation, GenesStates):    
    def __init__(
            self, path_to_tree, path_to_states, outdir, 
            gcode=2, db_mode="dict", path_to_db=None,
            rewrite_db=None, use_proba=False, proba_cutoff=0.05, 
            use_phylocoef=False, syn=False, syn_c=False, syn4f=False, derive_spectra=True,
            path_to_rates=None, cat_cutoff=2, save_exp_muts=False,
        ):
        for path in list(path_to_states) + [path_to_tree]:
            if not os.path.exists(path):
                raise ValueError(f"Path doesn't exist: {path}")

        CodonAnnotation.__init__(self, gcode)
        GenesStates.__init__(
            self, path_to_states, path_to_db, db_mode, rewrite_db, use_proba, 
            path_to_rates, cat_cutoff,
        )
        self.gcode = gcode
        self.use_proba = use_proba
        self.proba_cutoff = proba_cutoff
        self.use_phylocoef = use_phylocoef if use_proba else False
        self.derive_spectra = derive_spectra
        self.outdir = outdir
        self._save_exp_muts = save_exp_muts
        logger.info(f"Using gencode {gcode}")
        logger.info(f"Use probabilities of genomic states: {use_proba}")
        logger.info(f"Use phylogenetic uncertainty coefficient: {use_phylocoef}")
        logger.info(f"Derive spectra: {derive_spectra}")
        self.MUT_LABELS = ["all"]
        if syn:
            self.MUT_LABELS.append("syn")
        if syn_c:
            self.MUT_LABELS.append("syn_c")
        if syn4f:
            self.MUT_LABELS.append("ff")
        logger.info(f"Types of mutations to collect and process: {self.MUT_LABELS}")
        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        self.max_dist = self.fp_format(get_farthest_leaf(self.tree, 0.95))
        logger.info(
            f"Tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
            f"distance to farthest leaf: {self.max_dist: .2f}"
        )
        rnd_genome = self.get_random_genome()
        logger.info(f"Number of genes: {len(rnd_genome)}, number of sites: {[len(x) for x in rnd_genome.values()]}")
        if self.mask:
            logger.info(f"Number of invariable sites: {[len(x) - sum(x) for x in self.mask.values()]}")

    def open_handles(self, outdir):
        self.handle = dict()
        self.handle["mut"]  = open(os.path.join(outdir, "mutations.tsv"), "w")
        self.handle["freq"] = open(os.path.join(outdir, "expected_freqs.tsv"), "w")
        if self._save_exp_muts:
            self.handle["exp"] = open(os.path.join(outdir, "expected_mutations.tsv"), "w")

        if self.derive_spectra:
            self.handle["ms12"]   = open(os.path.join(outdir, "mutspec12.tsv"), "w")
            self.handle["ms192"]  = open(os.path.join(outdir, "mutspec192.tsv"), "w")
            self.handle["ms12s"]  = open(os.path.join(outdir, "mutspec12genes.tsv"), "w")
            self.handle["ms192g"] = open(os.path.join(outdir, "mutspec192genes.tsv"), "w")
        logger.info("Handles opened")

    def close_handles(self):
        for file in self.handle.values():
            file.close()
        logger.info("Handles closed")

    def extract_mutspec_from_tree(self):
        self.open_handles(self.outdir)

        logger.info("Start mutation extraction from tree")
        add_header = defaultdict(lambda: True)
        aln_size = self.genome_size
        dists_to_leafs = dict()
        visited_nodes = set()
        total_mut_num = 0
        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree), 1):
            if alt_node.name not in self.nodes:
                logger.warning(f"Skip edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{alt_node.name}' genome")
                continue
            if ref_node.name not in self.nodes:
                logger.warning(f"Skip edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{ref_node.name}' genome")
                continue

            # get genomes from storage
            ref_genome = self.get_genome(ref_node.name)
            alt_genome = self.get_genome(alt_node.name)

            # calculate phylogenetic uncertainty correction
            if self.use_phylocoef:
                if ref_node.name in dists_to_leafs:
                    dist_to_closest_leaf = dists_to_leafs[ref_node.name]
                else:
                    _, dist_to_closest_leaf = ref_node.get_closest_leaf()
                    dists_to_leafs[ref_node.name] = dist_to_closest_leaf

                dist_to_closest_leaf = self.fp_format(dist_to_closest_leaf)
                phylocoef = self.fp_format(1 - min(0.9999, dist_to_closest_leaf / self.max_dist))
            else:
                phylocoef = 1

            genome_nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_cxt_freqs  = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_mutations = []
            for gene in ref_genome:
                ref_seq = ref_genome[gene]
                alt_seq = alt_genome[gene]

                # collect state frequencies and
                # extract mutations and put in order columns
                if self.use_proba:
                    gene_exp_sbs12, gene_exp_sbs192 = self.collect_exp_mut_freqs_proba(ref_seq, phylocoef, self.mask[gene], self.MUT_LABELS, self.proba_cutoff)
                    gene_mut_df = self.extract_mutations_proba(ref_seq, alt_seq)
                    if self._save_exp_muts:
                        node_expected_sbs = self.collect_exp_muts_proba(ref_seq, phylocoef, self.mask[gene], self.MUT_LABELS, self.proba_cutoff)
                        node_expected_sbs["Node"] = ref_node.name
                        node_expected_sbs["Gene"] = gene
                else:
                    gene_exp_sbs12, gene_exp_sbs192 = self.collect_exp_mut_freqs(ref_seq, self.mask[gene], self.MUT_LABELS)
                    gene_mut_df = self.extract_mutations_simple(ref_seq, alt_seq)
                    if self._save_exp_muts:
                        node_expected_sbs = self.collect_exp_muts(ref_seq, self.mask[gene], self.MUT_LABELS)
                        node_expected_sbs["Node"] = ref_node.name
                        node_expected_sbs["Gene"] = gene

                # dump state frequencies
                if ref_node.name not in visited_nodes:
                    self.dump_expected_mutations(
                        gene_exp_sbs12, gene_exp_sbs192, ref_node.name, 
                        gene, self.handle["freq"], add_header["freqs"],
                    )
                    add_header["freqs"] = False

                    self.dump_table(node_expected_sbs, self.handle["exp"], add_header["exp"])
                    add_header["exp"] = False
                
                # summarize state frequencies over genome
                for lbl in self.MUT_LABELS:
                    for nucl, freq in gene_exp_sbs12[lbl].items():
                        genome_nucl_freqs[lbl][nucl] += freq
                    for trinucl, freq in gene_exp_sbs192[lbl].items():
                        genome_cxt_freqs[lbl][trinucl] += freq

                if gene_mut_df.shape[0] == 0:
                    continue
                
                if self.use_proba:
                    gene_mut_df["ProbaFull"] = phylocoef * gene_mut_df["ProbaMut"]

                gene_mut_df["RefNode"] = ref_node.name
                gene_mut_df["AltNode"] = alt_node.name
                gene_mut_df["Gene"] = gene
                
                # collect mutations of full genome
                genome_mutations.append(gene_mut_df)
                
                # calculate gene mutational spectra if there are at least 2 genes
                if self.derive_spectra and len(ref_genome) > 1:
                    for lbl in self.MUT_LABELS:
                        lbl_id = lbl2lbl_id("syn" if lbl == "syn_c" else lbl)

                        mutspec12 = calculate_mutspec(
                            gene_mut_df[gene_mut_df.Label >= lbl_id], gene_exp_sbs12[lbl], 
                            use_context=False, use_proba=self.use_proba
                        )
                        mutspec12["AltNode"] = alt_node.name
                        mutspec12["Label"] = lbl
                        mutspec12["Gene"]  = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec12, self.handle["ms12s"], add_header["ms12g"])
                        add_header["ms12g"] = False

                        if len(gene_mut_df) > 100:
                            mutspec192 = calculate_mutspec(
                                gene_mut_df[gene_mut_df.Label >= lbl_id], gene_exp_sbs192[lbl], 
                                use_context=True, use_proba=self.use_proba
                            )
                            mutspec192["RefNode"] = ref_node.name
                            mutspec192["AltNode"] = alt_node.name
                            mutspec192["Label"] = lbl
                            mutspec192["Gene"] = gene
                            # Dump gene mutspecs 
                            self.dump_table(mutspec192, self.handle["ms192g"], add_header["ms192g"])
                            add_header["ms192g"] = False

            visited_nodes.add(ref_node.name)
            
            if len(genome_mutations) == 0:
                logger.info(f"0 mutations from {ei:03} branch ({ref_node.name} - {alt_node.name})")
                continue

            genome_mutations_df = pd.concat(genome_mutations)
            del genome_mutations
            
            mut_num = genome_mutations_df.ProbaFull.sum() if self.use_proba else len(genome_mutations_df)
            total_mut_num += mut_num
            logger.info(f"{mut_num:.3f} mutations from branch {ei:03} ({ref_node.name} - {alt_node.name})")
            if mut_num > aln_size * 0.1:
                logger.warning(f"Observed too many mutations ({mut_num} > {aln_size} * 0.1) for branch ({ref_node.name} - {alt_node.name})")

            # dump mutations
            self.dump_table(genome_mutations_df, self.handle["mut"], add_header["mut"])
            add_header["mut"] = False
            
            # calculate full genome mutational spectra for all labels
            if self.derive_spectra:
                for lbl in self.MUT_LABELS:
                    lbl_id = lbl2lbl_id("syn" if lbl == "syn_c" else lbl)

                    mutspec12 = calculate_mutspec(
                        genome_mutations_df[genome_mutations_df.Label >= lbl_id],
                        genome_nucl_freqs[lbl], use_context=False, use_proba=self.use_proba
                    )
                    mutspec12["RefNode"] = ref_node.name
                    mutspec12["AltNode"] = alt_node.name
                    mutspec12["Label"] = lbl
                    mutspec192 = calculate_mutspec(
                        genome_mutations_df[genome_mutations_df.Label >= lbl_id],
                        genome_cxt_freqs[lbl], use_context=True, use_proba=self.use_proba
                    )
                    mutspec192["RefNode"] = ref_node.name
                    mutspec192["AltNode"] = alt_node.name
                    mutspec192["Label"] = lbl

                    # Dump genome spectra
                    self.dump_table(mutspec12,  self.handle["ms12"],  add_header["ms"])
                    self.dump_table(mutspec192, self.handle["ms192"], add_header["ms"])
                    add_header["ms"] = False

        logger.info(f"Processed {ei} tree edges")
        logger.info(f"Observed {total_mut_num:.3f} substitutions")
        logger.info("Extraction of mutations from phylogenetic tree completed succesfully")
        
        self.close_handles()

    def extract_mutations_proba(self, g1: np.ndarray, g2: np.ndarray):
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
            for cdn1, mut_cxt1, proba1 in self.sample_context(pos, pic, g1, self.proba_cutoff):
                cdn1_str = "".join(cdn1)
                for cdn2, mut_cxt2, proba2 in self.sample_context(pos, pic, g2, self.proba_cutoff):
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

    def collect_exp_mut_freqs_proba(
            self, cds: np.ndarray, phylocoef: float, 
            mask: Iterable[Union[int, bool]] = None, 
            labels = ["all", "syn", "ff"],  proba_cutoff=0.05,
        ):
        n = len(cds)
        if mask is not None and len(mask) != n:
            msg = f"Mask (len = {len(mask)}) must have same lenght as cds (len = {n})"
            logger.error(msg)
            logger.info("Termination")
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
            for cdn_tuple, cxt, proba in self.sample_context(pos, pic, cds, proba_cutoff):
                cdn = "".join(cdn_tuple)
                proba *= phylocoef  # adjusted proba
                nuc = cxt[1]
                sbs12_pattern = nuc + ">" + "{}"
                sbs192_pattern = cxt[0] + "[" + nuc + ">{}]" + cxt[-1]

                if ("syn" in labels or "syn_c" in labels) and len(cdn) == 3:
                    syn_codons = self.get_syn_codons(cdn, pic)
                    if "syn" in labels:
                        for alt_cdn in syn_codons:
                            alt_nuc = alt_cdn[pic]
                            sbs12_freqs["syn"][sbs12_pattern.format(alt_nuc)] += proba
                            sbs192_freqs["syn"][sbs192_pattern.format(alt_nuc)] += proba
                    if "syn_c" in labels and len(syn_codons) > 0:
                        for alt_nuc in self.nucl_order:
                            sbs12_freqs["syn_c"][sbs12_pattern.format(alt_nuc)] += proba
                            sbs192_freqs["syn_c"][sbs192_pattern.format(alt_nuc)] += proba

                for alt_nuc in self.nucl_order:
                    if alt_nuc == nuc:
                        continue
                    if "all" in labels:
                        sbs12_freqs["all"][sbs12_pattern.format(alt_nuc)] += proba
                        sbs192_freqs["all"][sbs192_pattern.format(alt_nuc)] += proba
                    if "pos3" in labels and pic == 2:
                        sbs12_freqs["pos3"][sbs12_pattern.format(alt_nuc)] += proba
                        sbs192_freqs["pos3"][sbs192_pattern.format(alt_nuc)] += proba
                    if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                        sbs12_freqs["ff"][sbs12_pattern.format(alt_nuc)] += proba
                        sbs192_freqs["ff"][sbs192_pattern.format(alt_nuc)] += proba
                        
        return sbs12_freqs, sbs192_freqs

    def collect_exp_muts_proba(
            self, cds: np.ndarray, phylocoef: float, 
            mask: Iterable[Union[int, bool]] = None, 
            labels = ["all", "syn", "ff"],  proba_cutoff=0.05,
        ):
        n = len(cds)
        if mask is not None and len(mask) != n:
            msg = f"Mask (len = {len(mask)}) must have same lenght as cds (len = {n})"
            logger.error(msg)
            logger.info("Termination")
            raise ValueError(msg)

        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        assert 0 < phylocoef <= 1, "Evol coefficient must be between 0 and 1"

        labels = set(labels)
        data = []
        for pos in range(1, n - 1):
            if mask is not None and not mask[pos]:
                continue
            pic = pos % 3  # 0-based
            for cdn_tuple, cxt, proba in self.sample_context(pos, pic, cds, proba_cutoff):
                cdn = "".join(cdn_tuple)
                proba *= phylocoef  # adjusted proba
                nuc = cxt[1]
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
                                "Proba": proba,
                            })
                    if "syn_c" in labels and len(syn_codons) > 0:
                        for alt_nuc in self.nucl_order:
                            data.append({
                                "Pos": pos + 1, "Pic": pic + 1,
                                "Mut": sbs192_pattern.format(alt_nuc),
                                "Cdn": cdn, "Label": "syn_c",
                                "Proba": proba,
                            })

                for alt_nuc in self.nucl_order:
                    if alt_nuc == nuc:
                        continue
                    if "all" in labels:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "all",
                            "Proba": proba,
                        })
                    if "pos3" in labels and pic == 2:
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "pos3",
                            "Proba": proba,
                        })
                    if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                        data.append({
                            "Pos": pos + 1, "Pic": pic + 1,
                            "Mut": sbs192_pattern.format(alt_nuc),
                            "Cdn": cdn, "Label": "syn4f",
                            "Proba": proba,
                        })

        exp_sbs = pd.DataFrame(data)
        return exp_sbs

    def sample_context(self, pos, pic, genome: np.ndarray, cutoff=0.01):
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

    @staticmethod
    def dump_table(df: pd.DataFrame, handle, header=False):
        if header:
            handle.write("\t".join(df.columns) + "\n")
        handle.write(df.to_csv(sep="\t", index=None, header=None))
    
    def dump_expected_mutations(self, gene_exp_sbs12, gene_exp_sbs192, node, gene, handle, header=False):
        if header:
            sbs12  = "\t".join(possible_sbs12)
            sbs192 = "\t".join(possible_sbs192)
            header = f"Node\tGene\tLabel\t{sbs12}\t{sbs192}\n"
            handle.write(header)
        
        for lbl in self.MUT_LABELS:
            sbs12  = "\t".join([str(gene_exp_sbs12[lbl][x]) for x in possible_sbs12])
            sbs192 = "\t".join([str(gene_exp_sbs192[lbl][x]) for x in possible_sbs192])
            row = f"{node}\t{gene}\t{lbl}\t{sbs12}\t{sbs192}\n"
            handle.write(row)


@click.command("MutSpec calculator", help="TODO")
@click.option("--tree", "path_to_tree", required=True, type=click.Path(True), help="Path to phylogenetic tree to collect mutations from")
@click.option("--states", "path_to_states", required=True, multiple=True, type=click.Path(True), help="Path to states of each node in the tree. Could be passed several states files, example: '--states file1 --states file2'")
@click.option("--outdir", required=True, type=click.Path(exists=False, writable=True), help="Directory which will contain output files with mutations and mutspecs")
@click.option("--gencode", required=True, type=int, help="Genetic code number to use in mutations annotation. Use 2 for vertebrate mitochondrial genes")
@click.option("--syn",   is_flag=True, default=False, help="Process synonymous mutations (expectations will be calculated using possible syn mutations counts)")
@click.option("--syn_c", is_flag=True, default=False, help="Process synonymous mutations (expectations will be calculated using possible syn mutations context counts)")
@click.option("--syn4f", is_flag=True, default=False, help="Process synonymous fourfold mutations")
@click.option("--proba", is_flag=True, default=False, help="Use states probabilities while mutations collecting")
@click.option("--pcutoff", "proba_cutoff", default=0.05, show_default=True, type=float, help="Cutoff of tri/tetranucleotide state probability, states with lower values will not be used in mutation collecting")
@click.option("--phylocoef/--no-phylocoef", is_flag=True, default=True, show_default=True, help="Use or don't use phylogenetic uncertainty coefficient. Considered only with --proba")
@click.option("--no-mutspec", "no_spectra", is_flag=True, default=False, show_default=True, help="Don't calculate mutspec, only mutations extraction")
@click.option("--save-exp-muts", is_flag=True, default=False, show_default=True, help="Save possible (expected) mutations to table")
@click.option("--rates", "path_to_rates", default=None, type=click.Path(True), help="Path to rates from IQTREE2")
@click.option("--cat-cutoff", type=int, default=1, show_default=True, help="Minimal category in rates file considered as variable. Default value 1 indicates that sites with category less than 1 will not be used in expected mutations counts")
@click.option("--write_db", is_flag=True, help="Write sqlite3 database instead of using dictionary for states. Usefull if you have low RAM. Time expensive")
@click.option("--db_path", "path_to_db", type=click.Path(writable=True), default="/tmp/states.db", show_default=True, help="Path to database with states. Use only with --write_db")
@click.option("--rewrite_db",  is_flag=True, default=False, help="Rewrite existing states database. Use only with --write_db")  # drop argument, replace by question
@click.option('-f', '--force', is_flag=True, help="Rewrite existing output directory")
@click.option('-q', '--quiet', is_flag=True, help="Quiet mode, suppress printing to screen (stderr)")
@click.option("--config", default=None, type=click.Path(True), help="Path to log-config file")
def main(
        path_to_tree, path_to_states, outdir, 
        gencode, syn, syn_c ,syn4f, proba, proba_cutoff,
        write_db, path_to_db, rewrite_db, 
        phylocoef, no_spectra, save_exp_muts, 
        path_to_rates, cat_cutoff,
        force, quiet, config,
    ):

    if os.path.exists(outdir):
        if not force:
            answer = input(f"Delete existing directory '{outdir}'? [Y/n] ")
            if answer.upper() != "Y" and answer != "":
                print("Interapted", file=sys.stderr)
                exit(0)
        rmtree(outdir)
        print(f"Directory '{outdir}' deleted", file=sys.stderr)
    os.makedirs(outdir)

    global logger
    _log_lvl = "CRITICAL" if quiet else None
    logfile = os.path.join(outdir, "run.log")
    logger = load_logger(path=config, stream_level=_log_lvl, filename=logfile)

    logger.info(f"Writing logs to '{logfile}'")
    logger.debug("Command: " + " ".join(sys.argv))
    logger.debug(f"Output directory '{outdir}' created")

    db_mode = "db" if write_db else "dict"
    derive_spectra = not no_spectra
    
    MutSpec(
        path_to_tree, path_to_states, outdir, gcode=gencode, 
        db_mode=db_mode, path_to_db=path_to_db, rewrite_db=rewrite_db, 
        use_proba=proba, proba_cutoff=proba_cutoff, use_phylocoef=phylocoef,
        syn=syn, syn_c=syn_c, syn4f=syn4f, derive_spectra=derive_spectra, 
        path_to_rates=path_to_rates, cat_cutoff=cat_cutoff, 
        save_exp_muts=save_exp_muts,
    ).extract_mutspec_from_tree()

if __name__ == "__main__":
    main()
