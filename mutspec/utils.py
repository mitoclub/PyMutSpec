import cProfile
import pstats
import re
import os
import sys
import sqlite3
from collections import defaultdict
from queue import Queue
from typing import List, Set, Tuple, Union

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA
from ete3 import PhyloTree
import tqdm

from custom_logging import logger

#################
## CODON FUNC  ##
#################


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
        
    def get_mut_effect(self, codon1: str, codon2: str, pic: int):
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


#################
##  TREE FUNC  ##
#################

def node_parent(node):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


def iter_tree_edges(tree: PhyloTree):
    discovered_nodes = set()
    discovered_nodes.add(tree.name)
    Q = Queue()
    Q.put(tree)

    while not Q.empty():
        cur_node = Q.get()
        for child in cur_node.children:
            Q.put(child)

        if cur_node.name not in discovered_nodes:
            discovered_nodes.add(cur_node.name)
            alt_node = cur_node
            ref_node = node_parent(alt_node)
            yield ref_node, alt_node


def get_farthest_leaf(tree: PhyloTree, quantile=None):
    """
    TODO check if tree is rooted
    """
    if quantile is None:
        _, md = tree.get_farthest_leaf()
        return md
    elif isinstance(quantile, (float, int)) and 0 <= quantile <= 1:
        distances = []
        for leaf in tree.iter_leaves():
            d = tree.get_distance(leaf)
            distances.append(d)
        md = np.quantile(distances, quantile)
    else:
        raise TypeError(f"quantile must be int, float or None, got {type(quantile)}")
    return md


#################
## OTHER FUNC  ##
#################


class GenomeStates:
    """
    tips:
    - use mode="dict" if your mtDNA tree is small (~1500 nodes require ~350MB of RAM); 
    big trees with 100000 nodes require tens GB of RAM, therefore use mode="db"
    """
    def __init__(self, 
            path_to_anc=None, path_to_leaves=None, path_to_db=None, 
            mode="dict", rewrite=False, *args, **kwargs
        ):
        self.mode = mode
        logger.debug(f"Genome states storage mode = '{mode}'")
        self.path_to_db = path_to_db
        if mode == "dict":
            self._prepare_node2genome(path_to_anc, path_to_leaves)
        elif mode == "db":
            if path_to_db is None:
                raise ValueError("Pass the path to database to use or write it")
            else:
                self._prepare_db(path_to_anc, path_to_leaves, rewrite)
        else:
            raise ValueError("Mode must be 'dict' or 'db'")
    
    def get_genome(self, node: str):
        if self.mode == "dict":
            return self.node2genome[node]
        elif self.mode == "db":
            genome = defaultdict(list)
            cur = self.con.cursor()
            for row in cur.execute(f"""SELECT Part, p_A, p_C, p_G, p_T FROM states 
                WHERE Node='{node}'
                ORDER by Part, Site"""):
                part = row[0]
                state = row[1:]
                genome[part].append(state)
            genome = {k: np.array(v, dtype=np.float16) for k, v in genome.items()}
            return genome
        else:
            raise ValueError("mode must be 'dict' or 'db'")

    def _prepare_db(self, path_to_anc, path_to_leaves, rewrite=False):
        """sequentially read tsv and write to db"""
        done = False
        if os.path.exists(self.path_to_db):
            if rewrite:
                os.remove(self.path_to_db)
            else:
                logger.info(f"Loading existing database from {self.path_to_db}")
                done = True
        try:
            con = sqlite3.connect(self.path_to_db)
        except:
            raise ValueError("Cannot connect to database by given path")
        
        cur = con.cursor()
        self.con = con
        if done:
            self.nodes = self._get_nodes()
            return

        logger.info("Writing new database file")
        cur.execute('''CREATE TABLE states
               (Node TEXT, Part INTEGER, Site INTEGER, State TEXT, p_A REAL, p_C REAL, p_G REAL, p_T REAL)'''
        )
        for path in [path_to_leaves, path_to_anc]:
            if path is None:
                continue
            logger.info(f"Loading {path} ...")
            handle = open(path, "r")
            header = handle.readline().strip()
            if header != "Node\tPart\tSite\tState\tp_A\tp_C\tp_G\tp_T":
                handle.close()
                con.close()
                raise ValueError(f"Inappropriate type of table, expected another columns order,\ngot {repr(header)}")

            for line in tqdm.tqdm(handle, total=8652300):
                row = line.strip().split()
                query = "INSERT INTO states VALUES ('{}',{},{},'{}',{},{},{},{})".format(*row)
                cur.execute(query)

            con.commit()
            handle.close()

        self.nodes = self._get_nodes()
        # con.close()
    
    def _get_nodes(self):
        nodes = set()
        cur = self.con.cursor()
        for node in cur.execute("SELECT DISTINCT Node from states"):
            nodes.add(node[0])
        return nodes

    def _prepare_node2genome(self, path_to_anc, path_to_leaves, states_dtype=np.float16):
        dtypes = {
            "p_A":  states_dtype, "p_C": states_dtype, 
            "p_G":  states_dtype, "p_T": states_dtype,
            "Site": np.int32,     "Part": np.int8,
        }
        usecols = ["Node", "Part", "Site", "p_A", "p_C", "p_G", "p_T"]
        node2genome = defaultdict(dict)
        for path in [path_to_anc, path_to_leaves]:
            if path is None:
                continue
            logger.info(f"Loading {path}...")
            states = pd.read_csv(path, sep="\t", comment='#', usecols=usecols, dtype=dtypes)
            aln_sizes = states.groupby("Node").apply(len)
            assert aln_sizes.nunique() == 1, "uncomplete leaves state table"
            states = states.sort_values(["Node", "Part", "Site"])
            gr = states.groupby(["Node", "Part"])
            for (node, part), gene_pos_ids in gr.groups.items():
                gene_df = states.loc[gene_pos_ids]
                gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
                node2genome[node][part] = gene_states

        self.node2genome = node2genome
        self.nodes = set(node2genome.keys())


def profiler(_func=None, *, nlines=10):
    def decorator_profiler(func):
        def wrapper_profiler(*args, **kwargs):
            profile = cProfile.Profile()
            profile.enable()
            value = func(*args, **kwargs)
            profile.disable()
            ps = pstats.Stats(profile, stream=sys.stderr)
            print(f"{'-' * 30}\nFunction: {func.__name__}\n{'-' * 30}", file=sys.stderr)
            ps.sort_stats('cumtime', 'calls')
            ps.print_stats(nlines)
            return value
        return wrapper_profiler

    if _func is None:
        return decorator_profiler
    else:
        return decorator_profiler(_func)


possible_sbs12 = {
    'A>C', 'A>G', 'A>T',
    'C>A', 'C>G', 'C>T',
    'G>A', 'G>C', 'G>T',
    'T>A', 'T>C', 'T>G'
}

possible_sbs192 = {
    "A[A>C]A", "A[A>C]C", "A[A>C]G", "A[A>C]T", "C[A>C]A", "C[A>C]C", "C[A>C]G", "C[A>C]T", 
    "G[A>C]A", "G[A>C]C", "G[A>C]G", "G[A>C]T", "T[A>C]A", "T[A>C]C", "T[A>C]G", "T[A>C]T", 
    "A[A>G]A", "A[A>G]C", "A[A>G]G", "A[A>G]T", "C[A>G]A", "C[A>G]C", "C[A>G]G", "C[A>G]T", 
    "G[A>G]A", "G[A>G]C", "G[A>G]G", "G[A>G]T", "T[A>G]A", "T[A>G]C", "T[A>G]G", "T[A>G]T", 
    "A[A>T]A", "A[A>T]C", "A[A>T]G", "A[A>T]T", "C[A>T]A", "C[A>T]C", "C[A>T]G", "C[A>T]T", 
    "G[A>T]A", "G[A>T]C", "G[A>T]G", "G[A>T]T", "T[A>T]A", "T[A>T]C", "T[A>T]G", "T[A>T]T", 
    "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", 
    "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", 
    "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", 
    "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", 
    "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", 
    "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", 
    "A[G>A]A", "A[G>A]C", "A[G>A]G", "A[G>A]T", "C[G>A]A", "C[G>A]C", "C[G>A]G", "C[G>A]T", 
    "G[G>A]A", "G[G>A]C", "G[G>A]G", "G[G>A]T", "T[G>A]A", "T[G>A]C", "T[G>A]G", "T[G>A]T", 
    "A[G>C]A", "A[G>C]C", "A[G>C]G", "A[G>C]T", "C[G>C]A", "C[G>C]C", "C[G>C]G", "C[G>C]T", 
    "G[G>C]A", "G[G>C]C", "G[G>C]G", "G[G>C]T", "T[G>C]A", "T[G>C]C", "T[G>C]G", "T[G>C]T", 
    "A[G>T]A", "A[G>T]C", "A[G>T]G", "A[G>T]T", "C[G>T]A", "C[G>T]C", "C[G>T]G", "C[G>T]T", 
    "G[G>T]A", "G[G>T]C", "G[G>T]G", "G[G>T]T", "T[G>T]A", "T[G>T]C", "T[G>T]G", "T[G>T]T", 
    "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", 
    "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", 
    "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", 
    "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", 
    "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", 
    "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T", 
}

possible_codons = {
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", 
    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", 
    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", 
    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", 
    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", 
    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", 
    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", 
    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT", 
}


if __name__ == "__main__":
    gcode = 2
    ca = CodonAnnotation(gcode)
    print(ca.is_syn_codons("ATC", "ATG"))
    print(read_start_stop_codons(2))
    print(ca.get_syn_number("ATG", 2))

    # path_to_anc = "./data/example_birds/anc_kg_states_birds.tsv"
    # path_to_leaves = "./data/example_birds/leaves_birds_states.tsv"
    # path_to_db = './data/states.db'

    path_to_states = "./data/states_sample.tsv"
    path_to_db = './data/states_sample.db'

    gs = GenomeStates(path_to_states, path_to_db=path_to_db, mode="db", rewrite=False)
    print(gs.get_genome("Node2"))

