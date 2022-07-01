"""
pic is position in codon

"""

import os
from collections import defaultdict
import sys

import click
import numpy as np
import pandas as pd
from ete3 import PhyloTree

from utils import (
    iter_tree_edges, profiler, calculate_mutspec, get_farthest_leaf,
    CodonAnnotation, GenomeStates, possible_codons
)
from mutspec.utils.logging import load_logger

logger = None


class MutSpec(CodonAnnotation, GenomeStates):    
    def __init__(
            self, path_to_tree, path_to_states, out_dir, 
            gcode=2, db_mode="dict", path_to_db="/tmp/states.db", 
            rewrite_db=None, proba_mode=False, proba_cutoff=0.01, 
            pastml_test=False, syn=True, syn4f=False,
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
        self.proba_cutoff = proba_cutoff
        self.pastml_test = pastml_test
        self.MUT_LABELS = ["all"]
        if syn:
            self.MUT_LABELS.append("syn")
        if syn4f:
            self.MUT_LABELS.append("ff")
        logger.info(f"Will calculate such mutspecs: {self.MUT_LABELS}")
        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        self.max_dist = self.fp_format(get_farthest_leaf(self.tree))
        logger.info(
            f"tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
            f"max distance to leaf: {self.max_dist: .2f}"
        )
        self.open_handles(out_dir)
        self.extract_mutspec_from_tree()
        self.close_handles()

    def open_handles(self, out_dir):
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

    def extract_mutspec_from_tree(self):
        logger.info("Start mutation extraction from tree")
        add_header = defaultdict(lambda: True)
        aln_size = self.genome_size
        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree), 1):
            if ref_node.name not in self.nodes or alt_node.name not in self.nodes:
                logger.debug(f"Pass edge '{ref_node.name}'-'{alt_node.name}' due to absence of genome")
                continue
            logger.debug(f"Extracting mutations from {ref_node.name} to {alt_node.name}")

            # get genomes from storage
            ref_genome = self.get_genome(ref_node.name)
            alt_genome = self.get_genome(alt_node.name)

            # calculate phylogenetic uncertainty correction
            _, dist_to_closest_leaf = ref_node.get_closest_leaf()
            dist_to_closest_leaf = self.fp_format(dist_to_closest_leaf)
            evol_speed_coef = self.fp_format(1 - min(1, dist_to_closest_leaf / self.max_dist))

            genome_nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_cxt_freqs  = {lbl: defaultdict(self.fp_format) for lbl in self.MUT_LABELS}
            genome_mutations = []
            for gene in ref_genome:
                ref_seq = ref_genome[gene]
                alt_seq = alt_genome[gene]

                # collect state frequencies
                if self.proba_mode:
                    gene_nucl_freqs, gene_cxt_freqs = self.collect_obs_mut_freqs_ptoba(ref_seq, evol_speed_coef)
                else:
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

                # extract mutations and put in order columns
                if self.proba_mode:
                    gene_mut_df = self.extract_mutations_proba(ref_seq, alt_seq)
                else:
                    gene_mut_df = self.extract_mutations_simple(ref_seq, alt_seq)

                if gene_mut_df.shape[0] == 0:
                    continue
                
                if self.proba_mode:
                    if self.pastml_test:
                        gene_mut_df["ProbaFull"] = gene_mut_df["ProbaMut"]
                    else:
                        gene_mut_df["ProbaFull"] = evol_speed_coef * gene_mut_df["ProbaMut"]
                gene_mut_df["RefNode"] = ref_node.name
                gene_mut_df["AltNode"] = alt_node.name
                gene_mut_df["Gene"] = gene
                
                # collect mutations of full genome
                genome_mutations.append(gene_mut_df)

                # calculate gene mutational spectra for all labels
                for lbl in self.MUT_LABELS:
                    mutspec12 = calculate_mutspec(gene_mut_df, gene_nucl_freqs[lbl], label=lbl, use_context=False, use_proba=self.proba_mode)
                    mutspec12["RefNode"] = ref_node.name
                    mutspec12["AltNode"] = alt_node.name
                    mutspec12["Label"] = lbl
                    mutspec12["Gene"]  = gene
                    # Dump gene mutspecs 
                    self.dump_table(mutspec12, self.handle["ms12s"], add_header["ms12g"])
                    add_header["ms12g"] = False

                if len(gene_mut_df) > 100:
                    for lbl in self.MUT_LABELS:
                        mutspec192 = calculate_mutspec(gene_mut_df, gene_cxt_freqs[lbl], label=lbl, use_context=True, use_proba=self.proba_mode)
                        mutspec192["RefNode"] = ref_node.name
                        mutspec192["AltNode"] = alt_node.name
                        mutspec192["Label"] = lbl
                        mutspec192["Gene"] = gene
                        # Dump gene mutspecs 
                        self.dump_table(mutspec192, self.handle["ms192g"], add_header["ms192g"])
                        add_header["ms192g"] = False

            if len(genome_mutations) == 0:
                logger.warning(f"Observed 0 mutations for branch ({ref_node.name} - {alt_node.name})")
                continue

            genome_mutations_df = pd.concat(genome_mutations)
            del genome_mutations
            if self.proba_mode and genome_mutations_df.ProbaFull.sum() > aln_size * 0.1:
                logger.warning(f"Observed too many mutations ({genome_mutations_df.ProbaFull.sum()} > {aln_size} * 0.1) for branch ({ref_node.name} - {alt_node.name})")
            if not self.proba_mode and len(genome_mutations_df) > aln_size * 0.1:
                logger.warning(f"Observed too many mutations ({len(genome_mutations_df)} > {aln_size} * 0.1) for branch ({ref_node.name} - {alt_node.name})")

            # dump mutations
            self.dump_table(genome_mutations_df, self.handle["mut"], add_header["mut"])
            add_header["mut"] = False
            
            # calculate full genome mutational spectra for all labels
            for lbl in self.MUT_LABELS:
                mutspec12 = calculate_mutspec(genome_mutations_df, genome_nucl_freqs[lbl], label=lbl, use_context=False, use_proba=self.proba_mode)
                mutspec12["RefNode"] = ref_node.name
                mutspec12["AltNode"] = alt_node.name
                mutspec12["Label"] = lbl
                mutspec192 = calculate_mutspec(genome_mutations_df, genome_cxt_freqs[lbl], label=lbl, use_context=True, use_proba=self.proba_mode)
                mutspec192["RefNode"] = ref_node.name
                mutspec192["AltNode"] = alt_node.name
                mutspec192["Label"] = lbl

                # Dump genome mutspecs
                self.dump_table(mutspec12,  self.handle["ms12"],  add_header["ms"])
                self.dump_table(mutspec192, self.handle["ms192"], add_header["ms"])
                add_header["ms"] = False

            if ei % 100 == 0:
                logger.info(f"Processed {ei} tree edges")

        logger.info(f"Processed {ei} tree edges")
        logger.info("MutSpec extraction done")

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

    def collect_obs_mut_freqs_ptoba(self, genome: np.ndarray, evol_speed_coef: float,  proba_cutoff=0.001):
        n = len(genome)
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        assert 0 <= evol_speed_coef <= 1, "Evol coefficient must be between 0 and 1"

        nucl_freqs = {lbl: defaultdict(self.fp_format) for lbl in ("all", "syn", "ff")}
        cxt_freqs = {lbl: defaultdict(self.fp_format) for lbl in ("all", "syn", "ff")}

        for pos in range(1, n - 1):
            pic = pos % 3  # 0-based
            for cdn, cxt, proba in self.sample_context(pos, pic, genome, proba_cutoff):
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
    

@click.command("MutSpec calculator", help="")
@click.option("--tree", "path_to_tree", required=True, type=click.Path(True), help="Path to phylogenetic tree to collect mutations from")
@click.option("--states", "path_to_states", required=True, multiple=True, type=click.Path(True), help="Path to states of each node in the tree. Could be passed several states files, example: '--states file1 --states file2'")
@click.option("--outdir", required=True, type=click.Path(exists=False, writable=True), help="Directory which will contain output files with mutations and mutspecs")
@click.option("--gencode", required=True, type=int, help="Genetic code number to use in mutations annotation. Use 2 for vertebrate mitochondrial genes")
@click.option("--syn", is_flag=True, required=False, default=False, help="Calculate synonymous mutspec")
@click.option("--syn4f", is_flag=True, required=False, default=False, help="Calculate synonymous fourfold mutspec")
@click.option("--proba", is_flag=True, required=False, default=False, help="Use states probabilities while mutations collecting")
@click.option("--pcutoff", "proba_cutoff", required=False, default=0.01, show_default=True, type=float, help="Cutoff of tri/tetranucleotide state probability, states with lower values will not be used in mutation collecting")
@click.option("--write_db", required=False, type=click.Choice(['dict', 'db'], case_sensitive=False), show_default=True, default="dict", help="Write sqlite3 database instead of using dictionary for states. Usefull if you have low RAM. Time expensive")
@click.option("--db_path", "path_to_db", required=False, type=click.Path(writable=True), default="/tmp/states.db", show_default=True, help="Path to database with states. Used only with --write_db")
@click.option("--rewrite_db", is_flag=True, required=False, default=False, help="Rewrite existing states database. Used only with --write_db")
@click.option("--pastml", is_flag=True, required=False, default=False, help="Run probability approach without phylogenetic uncertainty coefficient. Used for pastml run")
@click.option('-v', '--verbose', "verbosity", count=True, help="Verbosity level = DEBUG")
def main(
        path_to_tree, path_to_states, outdir, 
        gencode, syn, syn4f, proba, proba_cutoff, 
        write_db, path_to_db, rewrite_db, 
        pastml, verbosity,
    ):
    global logger
    os.makedirs(outdir)
    _log_lvl = "DEBUG" if verbosity >= 1 else None
    logger = load_logger(stream_level=_log_lvl, filename=os.path.join(outdir, "run.log"))
    logger.info(f"Output directory '{outdir}' created")
    logger.debug("Command: " + " ".join(sys.argv))
    MutSpec(
        path_to_tree, path_to_states, outdir, gcode=gencode, 
        db_mode=write_db, path_to_db=path_to_db, rewrite_db=rewrite_db, 
        proba_mode=proba, proba_cutoff=proba_cutoff, pastml_test=pastml,
        syn=syn, syn4f=syn4f,
    )


if __name__ == "__main__":
    main()


# if __name__ == "__main__":
#     path_to_tree =   "./data/example_nematoda/anc.treefile.rooted"
#     path_to_leaves = "./data/example_nematoda/leaves_states_nematoda.tsv"
#     path_to_anc = "./data/example_nematoda/nematoda_anc_mf/genes_states.tsv"
#     out_dir = "./data/processed/nematoda"

#     out_dir_lbl=""
#     out_dir_lbl = datetime.now().strftime("%d-%m-%y-%H-%M-%S") + "_" + out_dir_lbl
#     out_dir = os.path.join(out_dir, out_dir_lbl)
    
#     main(path_to_tree, path_to_anc, path_to_leaves, "anc_mf")
