#!/usr/bin/env python3
"""
pic is position in codon

"""

import os
import sys
import time
import random
from shutil import rmtree
from typing import Dict
from warnings import warn
import multiprocessing as mp
from datetime import datetime
from collections import defaultdict

import numpy as np
import pandas as pd
from ete3 import PhyloTree
from pymutspec.annotation import (
    CodonAnnotation, iter_tree_edges, get_tree_len,
)
from pymutspec.utils import load_logger, basic_logger

logger = basic_logger()

nucl_order5 = ['A', 'C', 'G', 'T', '-']


def calc_phylocoefs(tree: PhyloTree):
    tree_len = get_tree_len(tree, 'geom_mean')
    logger.info(f'Tree len = {tree_len:.3f}')
    phylocoefs = {tree.name: 1 - min(0.999, tree.get_closest_leaf()[1] / tree_len)}
    for node in tree.iter_descendants():
        _closest, d = node.get_closest_leaf()
        phylocoefs[node.name] = 1 - min(0.99999, d / tree_len)
    logger.info(f'Phylocoefs range: [{min(phylocoefs.values()):.3f}, {max(phylocoefs.values()):.3f}]')
    return phylocoefs


class GenomeStates:
    def __init__(self, path_to_states, path_to_gappy_sites):
        for path in list([path_to_states, path_to_gappy_sites]):
            if not os.path.exists(path):
                raise ValueError(f"Path to states doesn't exist: '{path}'")

        logger.info(f'Reading states "{path_to_states}"')
        states = pd.read_parquet(path_to_states)

        nodes_arr = states['Node'].unique()
        self.node2id = dict(zip(nodes_arr, range(len(nodes_arr))))
        self.nodes = set(nodes_arr)
        logger.debug('Node2id prepared')

        logger.info(f'Loaded {len(self.nodes)} genomes')

        states['NodeId'] = states['Node'].map(self.node2id).astype(np.int32)
        states.drop('Node', axis=1, inplace=True)
        logger.debug('Node column deleted')
        states.set_index('NodeId', inplace=True)
        logger.debug('States indexed')

        self.genome_size = len(states.loc[0])
        logger.info(f'Total alignement size = {self.genome_size}')

        self.nongappy_sites = self.read_nongappy_sites(path_to_gappy_sites)
        self.states = states
        logger.info('States initialiaztion is done')

    def __contains__(self, item: str):
        return item in self.nodes

    def __getitem__(self, item: str):
        return self.get_genome(item)

    def read_nongappy_sites(self, path):
        logger.debug(f'Reading gappy sites "{path}"')
        gappy_sites = set(pd.read_csv(path, header=None).values.flatten().tolist())
        nongappy_sites = np.array([i for i in range(self.genome_size) \
                                if i not in gappy_sites]) + 1 # 1-based
        logger.debug(f'Number of used (non-gappy) sites = {len(nongappy_sites)}')
        return nongappy_sites

    def get_genome(self, node: str) -> pd.DataFrame:
        if node not in self.nodes:
            raise ValueError(f'Passed node "{node}" have no genome')
        node_id = self.node2id[node]
        genome = self.states.loc[node_id].set_index('Site').loc[self.nongappy_sites]
        return genome


class GenomeStatesTotal():
    def __init__(self, path_to_states1: str, path_to_states2: str, path_to_gappy_sites: str):
        gs1 = GenomeStates(path_to_states1, path_to_gappy_sites)
        gs2 = GenomeStates(path_to_states2, path_to_gappy_sites)
        
        assert gs1.genome_size == gs2.genome_size
        self.genome_size = gs1.genome_size
        
        self.nodes = gs1.nodes.union(gs2.nodes)
        self.gs1 = gs1
        self.gs2 = gs2

    def __contains__(self, item):
        return item in self.nodes

    def __getitem__(self, item):
        if item in self.gs1:
            genome = self.gs1[item]
        elif item in self.gs2:
            genome = self.gs2[item]
        else:
            raise ValueError('node does not exist')
        return genome

    def get_genome(self, node: str) -> pd.DataFrame:
        if node in self.gs1:
            genome = self.gs1[node]
        elif node in self.gs2:
            genome = self.gs2[node]
        else:
            raise ValueError('node does not exist')
        return genome

    def get_random_genome(self) -> pd.DataFrame:
        node = random.choice(list(self.nodes))
        return self[node]


class Branch:
    def __init__(self, index: int, ref_name: str, alt_name: str, 
                 ref_genome: pd.DataFrame, alt_genome: pd.DataFrame, 
                 phylocoef: float, proba_cutoff: float, genome_size: int, 
                 site2index: Dict[int, int]) -> None:
        self.index = index
        self.ref_name = ref_name
        self.alt_name = alt_name
        self.ref_genome = ref_genome
        self.alt_genome = alt_genome
        self.phylocoef = phylocoef
        self.proba_cutoff = proba_cutoff
        self.genome_size = genome_size
        self.site2index = site2index

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
            logger.warning(f"Neighbouring positions of substitution are different: {cxt1}, {cxt2}")
            ok = False
        return ok

    def process_branch(self):
        edge_mutations = self.extract_mutations_proba(
            self.ref_genome, self.alt_genome, self.phylocoef, self.proba_cutoff)
        
        edge_mutations["RefNode"] = self.ref_name
        edge_mutations["AltNode"] = self.alt_name

        mut_num = edge_mutations['ProbaFull'].sum() if 'ProbaFull' in edge_mutations.columns else len(edge_mutations)
        if mut_num == 0:
            logger.info(f"No mutations from branch {self.index:03} ({self.ref_name} - {self.alt_name})")
        elif mut_num > self.genome_size * 0.1:
            logger.warning(f"Observed too many mutations ({mut_num} > {self.genome_size} * 0.1) "
                           f"for branch ({self.ref_name} - {self.alt_name})")
        else:
            logger.info(f"{mut_num:.2f} mutations from branch {self.index:03} ({self.ref_name} - {self.alt_name})")

        return edge_mutations

    # @profile
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
            logger.warning(f'Cannot collect local contexts for position {pos}')
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


class MutSpec(CodonAnnotation):
    def __init__(
            self, path_to_tree, outdir, gcode=2,
            use_proba=False, proba_cutoff=0.25, use_phylocoef=False,
            genomes: GenomeStatesTotal = None, num_processes=8.
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
        logger.info(f"Using gencode {gcode}")
        logger.info(f"Use probabilities of genomic states: {use_proba}")
        logger.info(f"Use phylogenetic uncertainty coefficient: {use_phylocoef}")
        logger.info(f"Minimal probability for mutations to use: {proba_cutoff}")

        self.fp_format = np.float32
        self.tree = PhyloTree(path_to_tree, format=1)
        logger.info(
            f"Tree loaded, number of leaf nodes: {len(self.tree)}, "
            f"total number of nodes: {len(self.tree.get_cached_content())}, "
        )
        rnd_genome = self.genomes.get_random_genome()
        logger.info(f"Number of MSA sites: {len(rnd_genome)}")

         # mapping of MSA sites to truncated states ids
        self.site2index = dict(zip(rnd_genome.index, range(len(rnd_genome))))

    def open_handles(self, outdir):
        self.handle = dict()
        self.handle["mut"]  = open(os.path.join(outdir, "mutations.tsv"), "w")
        logger.debug("Handles opened")

    def close_handles(self):
        for file in self.handle.values():
            file.close()
        logger.debug("Handles closed")

    def extract_mutspec_from_tree(self):
        t = self.num_processes
        if not isinstance(t, int) or t < 1:
            raise ValueError('num_processes must be positive integer')

        if t == 1:
            self._derive_mutspec()
        elif t > 1:
            self._derive_mutspec_parallel()

    def _derive_mutspec(self):
        logger.info("Start mutation extraction from tree")
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
        logger.info(f"Processed {edge_data[1]} tree edges")
        logger.info(f"Observed {total_mut_num:.3f} substitutions")
        logger.info("Extraction of mutations from phylogenetic tree completed succesfully")

    def _derive_mutspec_parallel(self):
        logger.info("Start parallel mutation extraction from tree")

        with mp.Pool(processes=self.num_processes) as pool:
            genome_mutations_lst = pool.map(Branch.process_branch, self.iter_branches())

        genome_mutations = pd.concat(genome_mutations_lst)

        self.open_handles(self.outdir)
        self.dump_table(genome_mutations, self.handle["mut"], True)
        self.close_handles()

        total_mut_num = genome_mutations['ProbaFull'].sum() if self.use_proba and \
            'ProbaFull' in genome_mutations.columns else len(genome_mutations)
        logger.info(f"Observed {total_mut_num:.3f} substitutions")
        logger.info("Extraction of mutations from phylogenetic tree completed succesfully")

    # def iter_chunks(self, n):
    #     lst = self.get_edges_data()
    #     for i in range(0, len(lst), n):
    #         yield lst[i:i + n]

    def iter_branches(self):
        """yield (self, edge_id, ref_node_name, alt_node_name, ref_genome, alt_genome, phylocoef)"""
        # calculate phylogenetic uncertainty correction
        if self.use_phylocoef:
            phylocoefs = calc_phylocoefs(self.tree)

        for ei, (ref_node, alt_node) in enumerate(iter_tree_edges(self.tree), 1):
            if alt_node.name not in self.genomes.nodes:
                logger.warning(f"Skip edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{alt_node.name}' genome")
                continue
            if ref_node.name not in self.genomes.nodes:
                logger.warning(f"Skip edge '{ref_node.name}'-'{alt_node.name}' due to absence of '{ref_node.name}' genome")
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
            logger.debug(f'Genomes retrieved for branch {ei}')

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


def main():
    debug = False
    current_datetime = str(datetime.now()).split('.')[0].replace(' ', '_')
    if debug:
        path_to_tree = 'data/tiny_tree2.nwk'
        path_to_states = 'data/tiny_states2.parquet'
        path_to_alignment = 'data/tiny_alignment2.parquet'
        outdir = f'data/tmp_{current_datetime}'
    else:
        path_to_tree = 'data/ancrec.root.nwk'
        path_to_states = 'data/states2.parquet'
        path_to_alignment = 'data/alignment2.parquet'
        outdir = 'data/multicore_run_v4'

    path_to_gappy_sites = './data/gappy_sites.csv'
    gencode = 2
    proba = phylocoef = True
    proba_cutoff = 0.25
    force = False
    num_processes = 64
    
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
    logfile = os.path.join(outdir, "run.log")
    logger = load_logger(filename=logfile, stream_level='INFO')

    logger.info(f"Writing logs to '{logfile}'")
    logger.debug("Command: " + " ".join(sys.argv))
    logger.debug(f"Output directory '{outdir}' created")
    if debug:
        logger.warning('DEBUG MODE ACTIVATED')

    start_time = time.time()

    genomes = GenomeStatesTotal(path_to_states, path_to_alignment, path_to_gappy_sites)

    MutSpec(
        path_to_tree, outdir, gcode=gencode, 
        use_proba=proba, proba_cutoff=proba_cutoff, use_phylocoef=phylocoef,
        genomes=genomes, num_processes=num_processes,
    ).extract_mutspec_from_tree()

    end_time = time.time()
    execution_time = end_time - start_time
    logger.info(f"Script completed taking {execution_time:.2f} seconds ({execution_time/3600:.2f} hours)")


if __name__ == "__main__":
    main()
