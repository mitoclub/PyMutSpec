"""
pic is position in codon

"""

import os
from collections import defaultdict
import sys
from shutil import rmtree

import click
import numpy as np
import pandas as pd
from ete3 import PhyloTree

from mutspec_utils.annotation import CodonAnnotation, get_farthest_leaf, iter_tree_edges, calculate_mutspec, lbl2lbl_id
from mutspec_utils.io import GenomeStates
from mutspec_utils.constants import possible_sbs12, possible_sbs192
from mutspec_utils.utils import load_logger, profiler

logger = None


class MutSpec(CodonAnnotation, GenomeStates):    
    def __init__(
            self, path_to_tree, path_to_states, out_dir, 
            gcode=2, db_mode="dict", path_to_db="/tmp/states.db", 
            rewrite_db=None, proba_mode=False, proba_cutoff=0.01, 
            use_phylocoef=False, syn=False, syn4f=False, no_mutspec=False,
        ):
        for path in list(path_to_states) + [path_to_tree]:
            if not os.path.exists(path):
                raise ValueError(f"Path doesn't exist: {path}")

        CodonAnnotation.__init__(self, gcode)
        GenomeStates.__init__(
            self, path_to_states, path_to_db, db_mode, rewrite_db, proba_mode
        )
        self.gcode = gcode
        logger.info(f"Using gencode {gcode}")
        self.proba_mode = proba_mode
        logger.info(f"Use probabilities of genomic states: {proba_mode}")
        self.proba_cutoff = proba_cutoff
        self.use_phylocoef = use_phylocoef
        self.no_mutspec = no_mutspec
        if proba_mode:
            logger.info(f"Use phylogenetic uncertainty coefficient: {use_phylocoef}")
        self.MUT_LABELS = ["all"]
        if syn:
            self.MUT_LABELS.append("syn")
        if syn4f:
            self.MUT_LABELS.append("ff")
        logger.info(f"Types of mutations to collect and process: {self.MUT_LABELS}")
        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        self.max_dist = self.fp_format(get_farthest_leaf(self.tree))
        logger.info(
            f"Tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
            f"distance to farthest leaf: {self.max_dist: .2f}"
        )
        self.open_handles(out_dir)
        self.extract_mutspec_from_tree()
        self.close_handles()

    def open_handles(self, out_dir):
        self.handle = dict()
        path_to_mutations  = os.path.join(out_dir, "mutations.tsv")
        path_to_nucl_freqs = os.path.join(out_dir, "expected_mutations.tsv")
        self.handle["mut"] = open(path_to_mutations, "w")
        self.handle["freq"]   = open(path_to_nucl_freqs, "w")
        if not self.no_mutspec:
            path_to_mutspec12  = os.path.join(out_dir, "mutspec12.tsv")
            path_to_mutspec192 = os.path.join(out_dir, "mutspec192.tsv")
            path_to_mutspec_genes12  = os.path.join(out_dir, "mutspec12genes.tsv")
            path_to_mutspec_genes192 = os.path.join(out_dir, "mutspec192genes.tsv")
            self.handle["ms12"]   = open(path_to_mutspec12, "w")
            self.handle["ms192"]  = open(path_to_mutspec192, "w")
            self.handle["ms12s"]  = open(path_to_mutspec_genes12, "w")
            self.handle["ms192g"] = open(path_to_mutspec_genes192, "w")
        logger.info("Handles opened")

    def close_handles(self):
        for file in self.handle.values():
            file.close()
        logger.info("Handles closed")

    def extract_mutspec_from_tree(self):
        logger.info("Start mutation extraction from tree")
        add_header = defaultdict(lambda: True)
        aln_size = self.genome_size
        dists_to_leafs = dict()
        total_mut_num = 0
        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree), 1):
            if alt_node.name not in self.nodes:
                logger.warning(f"Pass edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{alt_node.name}' genome")
                continue
            if ref_node.name not in self.nodes:
                logger.warning(f"Pass edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{ref_node.name}' genome")
                continue

            # get genomes from storage
            ref_genome = self.get_genome(ref_node.name)
            alt_genome = self.get_genome(alt_node.name)

            # calculate phylogenetic uncertainty correction
            if ref_node.name in dists_to_leafs:
                dist_to_closest_leaf = dists_to_leafs[ref_node.name]
            else:
                _, dist_to_closest_leaf = ref_node.get_closest_leaf()

            dist_to_closest_leaf = self.fp_format(dist_to_closest_leaf)
            phylocoef = self.fp_format(1 - min(1, dist_to_closest_leaf / self.max_dist))

            genome_nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_cxt_freqs  = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_mutations = []
            for gene in ref_genome:
                ref_seq = ref_genome[gene]
                alt_seq = alt_genome[gene]

                # collect state frequencies
                if self.proba_mode:
                    gene_exp_sbs12, gene_exp_sbs192 = self.collect_exp_mut_freqs_proba(ref_seq, phylocoef, self.MUT_LABELS)
                else:
                    gene_exp_sbs12, gene_exp_sbs192 = self.collect_exp_mut_freqs(ref_seq)

                # dump state frequencies
                if ref_node.name not in dists_to_leafs:
                    self.dump_expected_mutations(
                        gene_exp_sbs12, gene_exp_sbs192, ref_node.name, 
                        gene, self.handle["freq"], add_header["freqs"],
                    )
                    add_header["freqs"] = False
                
                # summarize state frequencies over genome
                for lbl in self.MUT_LABELS:
                    for nucl, freq in gene_exp_sbs12[lbl].items():
                        genome_nucl_freqs[lbl][nucl] += freq
                    for trinucl, freq in gene_exp_sbs192[lbl].items():
                        genome_cxt_freqs[lbl][trinucl] += freq

                # extract mutations and put in order columns
                if self.proba_mode:
                    gene_mut_df = self.extract_mutations_proba(ref_seq, alt_seq)
                else:
                    gene_mut_df = self.extract_mutations_simple(ref_seq, alt_seq)

                if gene_mut_df.shape[0] == 0:
                    continue
                
                if self.proba_mode:
                    if self.use_phylocoef:
                        gene_mut_df["ProbaFull"] = phylocoef * gene_mut_df["ProbaMut"]
                    else:
                        gene_mut_df["ProbaFull"] = gene_mut_df["ProbaMut"]
                gene_mut_df["RefNode"] = ref_node.name
                gene_mut_df["AltNode"] = alt_node.name
                gene_mut_df["Gene"] = gene
                
                # collect mutations of full genome
                genome_mutations.append(gene_mut_df)
                
                # calculate gene mutational spectra for all labels
                if not self.no_mutspec:
                    for lbl in self.MUT_LABELS:
                        mutspec12 = calculate_mutspec(
                            gene_mut_df[gene_mut_df.Label >= lbl2lbl_id(lbl)], gene_exp_sbs12[lbl], 
                            use_context=False, use_proba=self.proba_mode
                        )
                        mutspec12["AltNode"] = alt_node.name
                        mutspec12["Label"] = lbl
                        mutspec12["Gene"]  = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec12, self.handle["ms12s"], add_header["ms12g"])
                        add_header["ms12g"] = False

                        if len(gene_mut_df) > 100:
                            mutspec192 = calculate_mutspec(
                                gene_mut_df[gene_mut_df.Label >= lbl2lbl_id(lbl)], gene_exp_sbs192[lbl], 
                                use_context=True, use_proba=self.proba_mode
                            )
                            mutspec192["RefNode"] = ref_node.name
                            mutspec192["AltNode"] = alt_node.name
                            mutspec192["Label"] = lbl
                            mutspec192["Gene"] = gene
                            # Dump gene mutspecs 
                            self.dump_table(mutspec192, self.handle["ms192g"], add_header["ms192g"])
                            add_header["ms192g"] = False

            # save current node distance to closest leaf after genes processing
            if ref_node.name not in dists_to_leafs:
                dists_to_leafs[ref_node.name] = dist_to_closest_leaf
            
            if len(genome_mutations) == 0:
                logger.info(f"Observed 0 mutations for branch ({ref_node.name} - {alt_node.name})")
                continue

            genome_mutations_df = pd.concat(genome_mutations)
            del genome_mutations
            
            mut_num = genome_mutations_df.ProbaFull.sum() if self.proba_mode else len(genome_mutations_df)
            total_mut_num += mut_num
            logger.info(f"Observed {mut_num:.3f} mutations for branch ({ref_node.name} - {alt_node.name})")
            if mut_num > aln_size * 0.1:
                logger.warning(f"Observed too many mutations ({mut_num} > {aln_size} * 0.1) for branch ({ref_node.name} - {alt_node.name})")

            # dump mutations
            self.dump_table(genome_mutations_df, self.handle["mut"], add_header["mut"])
            add_header["mut"] = False
            
            # calculate full genome mutational spectra for all labels
            if not self.no_mutspec:
                for lbl in self.MUT_LABELS:
                    mutspec12 = calculate_mutspec(
                        genome_mutations_df[genome_mutations_df.Label >= lbl2lbl_id(lbl)],
                        genome_nucl_freqs[lbl], use_context=False, use_proba=self.proba_mode
                    )
                    mutspec12["RefNode"] = ref_node.name
                    mutspec12["AltNode"] = alt_node.name
                    mutspec12["Label"] = lbl
                    mutspec192 = calculate_mutspec(
                        genome_mutations_df[genome_mutations_df.Label >= lbl2lbl_id(lbl)],
                        genome_cxt_freqs[lbl], use_context=True, use_proba=self.proba_mode
                    )
                    mutspec192["RefNode"] = ref_node.name
                    mutspec192["AltNode"] = alt_node.name
                    mutspec192["Label"] = lbl

                    # Dump genome spectra
                    self.dump_table(mutspec12,  self.handle["ms12"],  add_header["ms"])
                    self.dump_table(mutspec192, self.handle["ms192"], add_header["ms"])
                    add_header["ms"] = False

            if ei % 100 == 0:
                logger.info(f"Processed {ei} tree edges")

        logger.info(f"Processed {ei} tree edges")
        logger.info(f"Observed {total_mut_num:.3f} substitutions")
        logger.info("Extraction of mutations from phylogenetic tree completed succesfully")

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
        - nucl_freqs - dict[lbl: dict[{ACGT}: int]] - nucleotide frequencies for all, syn and ff positions
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

    def collect_exp_mut_freqs_proba(self, genome: np.ndarray, phylocoef: float, labels = ["all", "syn", "ff"],  proba_cutoff=0.001):
        n = len(genome)
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        assert 0 <= phylocoef <= 1, "Evol coefficient must be between 0 and 1"

        sbs12_freqs = {lbl: defaultdict(int) for lbl in labels}
        sbs192_freqs = {lbl: defaultdict(int) for lbl in labels}
        labels = set(labels)

        for pos in range(1, n - 1):
            pic = pos % 3  # 0-based
            for cdn_tuple, cxt, proba in self.sample_context(pos, pic, genome, proba_cutoff):
                cdn = "".join(cdn_tuple)
                proba *= phylocoef  # adjusted proba
                nuc = cxt[1]
                mut_base12 = nuc + ">" + "{}"
                mut_base192 = cxt[0] + "[" + nuc + ">{}]" + cxt[-1]

                if "syn" in labels:
                    syn_codons = self.get_syn_codons(cdn, pic)
                    for alt_cdn in syn_codons:
                        alt_nuc = alt_cdn[pic]
                        sbs12_freqs["syn"][mut_base12.format(alt_nuc)] += proba
                        sbs192_freqs["syn"][mut_base192.format(alt_nuc)] += proba

                for alt_nuc in self.nucl_order:
                    if alt_nuc == nuc:
                        continue
                    if "all" in labels:
                        sbs12_freqs["all"][mut_base12.format(alt_nuc)] += proba
                        sbs192_freqs["all"][mut_base192.format(alt_nuc)] += proba
                    if "pos3" in labels and pic == 2:
                        sbs12_freqs["pos3"][mut_base12.format(alt_nuc)] += proba
                        sbs192_freqs["pos3"][mut_base192.format(alt_nuc)] += proba
                    if "ff" in labels and pic == 2 and self.is_fourfold(cdn):
                        sbs12_freqs["ff"][mut_base12.format(alt_nuc)] += proba
                        sbs192_freqs["ff"][mut_base192.format(alt_nuc)] += proba
                        
        return sbs12_freqs, sbs192_freqs

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
@click.option("--syn", is_flag=True, default=False, help="Calculate synonymous mutspec")
@click.option("--syn4f", is_flag=True, default=False, help="Calculate synonymous fourfold mutspec")
@click.option("--proba", is_flag=True, default=False, help="Use states probabilities while mutations collecting")
@click.option("--pcutoff", "proba_cutoff", default=0.01, show_default=True, type=float, help="Cutoff of tri/tetranucleotide state probability, states with lower values will not be used in mutation collecting")
@click.option("--phylocoef/--no-phylocoef", is_flag=True, default=True, show_default=True, help="Use or don't use phylogenetic uncertainty coefficient. Use only with --proba")
@click.option("--no-mutspec", is_flag=True, default=False, show_default=True, help="Don't calculate mutspec, only mutations extraction")
@click.option("--write_db", type=click.Choice(['dict', 'db'], case_sensitive=False), show_default=True, default="dict", help="Write sqlite3 database instead of using dictionary for states. Usefull if you have low RAM. Time expensive")
@click.option("--db_path", "path_to_db", type=click.Path(writable=True), default="/tmp/states.db", show_default=True, help="Path to database with states. Use only with --write_db")
@click.option("--rewrite_db",  is_flag=True, default=False, help="Rewrite existing states database. Use only with --write_db")
@click.option('-f', '--force', is_flag=True, help="Rewrite existing output directory")
@click.option('-v', '--verbose', "verbosity", count=True, help="Verbosity level = DEBUG")
def main(
        path_to_tree, path_to_states, outdir, 
        gencode, syn, syn4f, proba, proba_cutoff, 
        write_db, path_to_db, rewrite_db, 
        phylocoef, no_mutspec, force, verbosity,
    ):
    if os.path.exists(outdir):
        if not force:
            answer = input(f"Delete existing directory '{outdir}'? [Y/n] ")
            print(repr(answer), answer == "")
            if answer.upper() != "Y" and answer != "":
                print("Interapted")
                exit(0)
        rmtree(outdir)
        print(f"Existing output directory '{outdir}' deleted")
    os.makedirs(outdir)
    print(f"Output directory '{outdir}' created")
    global logger
    _log_lvl = "DEBUG" if verbosity >= 1 else None
    logfile = os.path.join(outdir, "run.log")
    logger = load_logger(stream_level=_log_lvl, filename=logfile)
    logger.info(f"Writing logs to '{logfile}'")
    logger.debug("Command: " + " ".join(sys.argv))
    MutSpec(
        path_to_tree, path_to_states, outdir, gcode=gencode, 
        db_mode=write_db, path_to_db=path_to_db, rewrite_db=rewrite_db, 
        proba_mode=proba, proba_cutoff=proba_cutoff, use_phylocoef=phylocoef,
        syn=syn, syn4f=syn4f, no_mutspec=no_mutspec,
    )


if __name__ == "__main__":
    main()
