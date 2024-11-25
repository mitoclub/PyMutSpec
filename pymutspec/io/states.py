import os
import sys
import random
import sqlite3
from collections import defaultdict
from typing import List

import tqdm
import numpy as np
import pandas as pd
from Bio import SeqIO

from ..utils import basic_logger


class GenomeStates:
    def __init__(self, path_to_states, path_to_gappy_sites=None, format='tsv', logger=None):
        for path in list([path_to_states, path_to_gappy_sites]):
            if not os.path.exists(path):
                raise ValueError(f"Path to states doesn't exist: '{path}'")

        self.logger = logger or basic_logger()

        self.logger.info(f'Reading states "{path_to_states}"')
        states = self.read_states(path, format)

        nodes_arr = states['Node'].unique()
        self.node2id = dict(zip(nodes_arr, range(len(nodes_arr))))
        self.nodes = set(nodes_arr)
        self.logger.info(f'Loaded {len(self.nodes)} genomes')

        states['NodeId'] = states['Node'].map(self.node2id)
        if len(self.nodes) < 2 ** 31 - 1:
            states['NodeId'] = states['NodeId'].astype(np.int32)

        states.drop('Node', axis=1, inplace=True)
        self.logger.debug('Node column deleted')

        states.set_index('NodeId', inplace=True)
        self.logger.debug('States reindexed')

        self.genome_size = len(states.loc[0])
        self.logger.info(f'Total alignement size = {self.genome_size}')

        self.nongappy_sites = self.read_nongappy_sites(path_to_gappy_sites)
        self.states = states
        self.logger.info('States initialiaztion is done')

    def __contains__(self, item: str):
        return item in self.nodes

    def __getitem__(self, item: str):
        return self.get_genome(item)

    def read_states(self, path, fmt='iqtree'):
        fpt = np.float32
        dtype = {
            "p_A":  fpt, "p_C": fpt, "p_G":  fpt, "p_T": fpt,
            "Site": np.int32, "Node": str,
        }
        if fmt == 'iqtree':
            states = pd.read_csv(path, sep='\t', comment='#', dtype=dtype)
        if fmt == 'csv':
            states = pd.read_csv(path, dtype=dtype)
        elif fmt == 'tsv':
            states = pd.read_csv(path, sep='\t', dtype=dtype)
        elif fmt == 'parquet':
            states = pd.read_parquet(path)
        return states

    def read_nongappy_sites(self, path):
        '''path must contain list of 1-based sites of alignment that must be excluded from loaded alignment
        
        return array of sorted non-gappy sites in 1-based mode
        '''
        if path is None:
            nongappy_sites = np.arange(1, self.genome_size + 1)
            self.logger.debug('Use all alignment sites')
        else:
            self.logger.debug(f'Reading gappy sites "{path}"')
            gappy_sites = set(pd.read_csv(path, header=None).values.flatten().tolist())
            nongappy_sites = np.array([i for i in range(1, self.genome_size + 1) \
                                    if i not in gappy_sites])
            self.logger.debug(f'Number of used (non-gappy) sites = {len(nongappy_sites)}')
        return nongappy_sites

    def get_genome(self, node: str) -> pd.DataFrame:
        if node not in self.nodes:
            raise ValueError(f'Input node "{node}" have no genome')
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


# class GenesStates:
#     """
#     Load genome (gene) states files and access gene by record name
    
#     tips:
#     - use mode="dict" if your mtDNA tree is small (~1500 nodes x 10,000 nucleotides require ~350MB of RAM); 
#     big trees with 100,000 nodes require tens of GB of RAM, therefore use mode="db"

#     Arguments
#     ---------
#     path_states: list or iterable of string
#         files with states; if some records not presented in all files table will be filled by '-' (gaps) 
#     TODO
#     states_fmt: string
#         format of files aligment; supported any format from biopython (fasta, phylip, etc.)
    
#     Return
#     ---------
#         states_obj: GenesStates - 
    
#     """
#     def __init__(
#             self, path_states: List[str], path_to_db=None, mode="dict", rewrite=False, 
#             use_proba=True, path_to_rates=None, cat_cutoff=1, states_fmt="table",
#         ):
#         fmt_variants = ["fasta", "phylip", "table"]
#         if states_fmt not in fmt_variants:
#             raise ValueError(f"Appropriate states_fmt: {repr(fmt_variants)}")
#         if states_fmt != "table" and use_proba == True:
#             print("Ignore use_proba option and forcely set use_proba=False due to input alignment that doesn't contain probabilities", file=sys.stderr)
#             use_proba = False
        
#         self.mode = mode
#         self.input_states_fmt = states_fmt
#         self.use_proba = use_proba
#         print(f"Genome states storage mode = '{mode}'", file=sys.stderr)
#         self.path_to_db = path_to_db
#         if mode == "dict":
#             if states_fmt == "table":
#                 self._prepare_node2genome(path_states)
#             else:
#                 states = self.read_alignment(path_states)
#                 self._prepare_node2genome(states)
#         elif mode == "db":
#             if states_fmt in ["fasta", "phylip"]:
#                 raise NotImplementedError
            
#             if path_to_db is None:
#                 raise ValueError("Pass the path to database to use or write it")
#             else:
#                 self._prepare_db(path_states, rewrite)
#         else:
#             raise ValueError("Mode must be 'dict' or 'db'")
        
#         genes_sizes = {part: len(x) for part, x in self.get_random_genome().items()}
#         self.genome_size = sum(genes_sizes.values())
#         self.category = self.read_rates(path_to_rates) if path_to_rates else None
#         self.mask = None

#         if self.category: 
#             for part in genes_sizes:
#                 if genes_sizes[part] != len(self.category[part]):
#                     raise RuntimeError("Wrong rates: number of positions in alignment are not equal to number of rates", file=sys.stderr)
#             self.mask = self.get_mask(self.category, cat_cutoff)
    
#     def get_genome(self, node: str):
#         if self.mode == "dict":
#             return self.node2genome[node]
            
#         elif self.mode == "db":
#             genome_raw = defaultdict(list)
#             cur = self.con.cursor()
#             if self.use_proba:
#                 query = f"""SELECT Part, Site, p_A, p_C, p_G, p_T FROM states WHERE Node='{node}'"""
#                 dtype = [
#                     ("Site", np.int32),
#                     ("p_A", np.float32), ("p_C", np.float32), 
#                     ("p_G", np.float32), ("p_T", np.float32),
#                 ]
#             else:
#                 query = f"""SELECT Part, Site, State FROM states WHERE Node='{node}'"""
#                 dtype = [("Site", np.int32), ("State", np.object_)]

#             for row in cur.execute(query):
#                 part = str(row[0])
#                 state = row[1:]
#                 genome_raw[part].append(state)
            
#             genome = {}
#             for part, state_tup in genome_raw.items():
#                 state = np.array(state_tup, dtype=dtype)
#                 state_sorted = np.sort(state, order="Site")
#                 state_sorted_devoid = np.array(list(map(list, state_sorted)))
#                 state = state_sorted_devoid[:, 1:] if self.use_proba else state_sorted_devoid[:, 1]
#                 genome[part] = state
#             return genome

#         else:
#             raise ValueError("mode must be 'dict' or 'db'")
    
#     def get_random_genome(self):
#         node = next(iter(self.nodes))
#         return self.get_genome(node)

#     def _prepare_db(self, path_states, rewrite=False):
#         """sequentially read tsv and write to db"""
#         done = False
#         if os.path.exists(self.path_to_db):
#             if rewrite:
#                 os.remove(self.path_to_db)
#             else:
#                 print(f"Loading existing database from {self.path_to_db}", file=sys.stderr)
#                 done = True
#         try:
#             con = sqlite3.connect(self.path_to_db)
#         except:
#             raise ValueError("Cannot connect to database by given path")
        
#         cur = con.cursor()
#         self.con = con
#         if done:
#             self.nodes = self._get_nodes_from_db()
#             return

#         print("Writing new database file", file=sys.stderr)
#         cur.execute('''CREATE TABLE states
#                (Node TEXT, Part INTEGER, Site INTEGER, State TEXT, p_A REAL, p_C REAL, p_G REAL, p_T REAL)'''
#         )
#         for path in path_states:
#             if path is None:
#                 continue
#             print(f"Loading {path} ...", file=sys.stderr)
#             handle = open(path, "r")
#             header = handle.readline().strip()
#             if header != "Node\tPart\tSite\tState\tp_A\tp_C\tp_G\tp_T":
#                 handle.close()
#                 con.close()
#                 raise ValueError(f"Inappropriate type of table, expected another columns order,\ngot {repr(header)}")

#             for line in tqdm.tqdm(handle, total=8652300):  # TODO estimate total
#                 row = line.strip().split()
#                 query = "INSERT INTO states VALUES ('{}',{},{},'{}',{},{},{},{})".format(*row)
#                 cur.execute(query)

#             con.commit()
#             handle.close()

#         self.nodes = self._get_nodes_from_db()
#         # con.close()
    
#     def _get_nodes_from_db(self):
#         nodes = set()
#         cur = self.con.cursor()
#         for node in cur.execute("SELECT DISTINCT Node from states"):
#             nodes.add(node[0])
#         return nodes
    
#     def states2dct(self, states: pd.DataFrame, out=None):
#         # all genes (genomes) must have same length
#         aln_sizes = states.groupby("Node").apply(len)
#         assert aln_sizes.nunique() == 1, "uncomplete state table: some site states absent in some node genomes"
        
#         node2genome = defaultdict(dict) if out is None else out
#         gr = states.sort_values(["Node", "Part", "Site"]).groupby(["Node", "Part"])
#         for (node, part), gene_df in gr:
#             if self.use_proba:
#                 gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
#             else:
#                 gene_states = gene_df.State.values
#             node2genome[node][part] = gene_states
#         return node2genome

#     def _prepare_node2genome(self, path_states, states_dtype=np.float32):
#         node2genome = defaultdict(dict)
#         if isinstance(path_states, pd.DataFrame):
#             states = path_states
#             self.states2dct(states, node2genome)
#         elif isinstance(path_states, (list, tuple)):
#             dtypes = {
#                 "p_A":  states_dtype, "p_C": states_dtype, 
#                 "p_G":  states_dtype, "p_T": states_dtype,
#                 "Site": np.int32, "Node": str, "Part": str,
#             }
#             if self.use_proba:
#                 usecols = ["Node", "Part", "Site", "p_A", "p_C", "p_G", "p_T"]
#             else:
#                 usecols = ["Node", "Part", "Site", "State"]

#             for path in path_states:
#                 if path is None:
#                     continue
#                 print(f"Loading {path}...", file=sys.stderr)
#                 states = pd.read_csv(path, sep="\t", comment='#', usecols=usecols, dtype=dtypes)
#                 self.states2dct(states, node2genome)

#                 # aln_sizes = states.groupby("Node").apply(len)
#                 # assert aln_sizes.nunique() == 1, "uncomplete leaves state table"
#                 # states = states.sort_values(["Node", "Part", "Site"])
#                 # gr = states.groupby(["Node", "Part"])
#                 # for (node, part), gene_pos_ids in tqdm.tqdm(gr.groups.items()):
#                 #     gene_df = states.loc[gene_pos_ids]  # TODO simplify: iterate over gr
#                 #     if self.use_proba:
#                 #         gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
#                 #     else:
#                 #         gene_states = gene_df.State.values
#                 #     node2genome[node][part] = gene_states
#         else:

#             raise NotImplementedError('path_to_states must be list of paths or pd.DataFrame')
        
#         self.node2genome = node2genome
#         self.nodes = set(node2genome.keys())

#     @staticmethod
#     def read_alignment(files: list, fmt="fasta"):
#         """
#         Read files alignments and prepare states table

#         Arguments
#         ---------
#         files: list or iterable of string
#             files with alignments; if some records not presented in all files table will be filled by '-' (gaps) 
#         fmt: string
#             format of files aligment; supported any format from biopython (fasta, phylip, etc.)
        
#         Return
#         ---------
#             states: pd.DataFrame
#         """
#         ngenes = len(files)  
#         columns = "Node Part Site State p_A p_C p_G p_T".split()
#         nucls = "ACGT"
#         aln_lens = dict()
#         files = set(files)
#         history = defaultdict(list)
#         data = []
#         for filepath in files:
#             part = "1" if ngenes == 1 else os.path.basename(filepath).replace(".fna", "")
#             alignment = SeqIO.parse(filepath, fmt)
#             for rec in alignment:
#                 node = rec.name
#                 history[node].append(part)
#                 seq = str(rec.seq)
#                 for site, state in enumerate(seq, 1):
#                     site_data = [node, part, site, state]
#                     for nucl in nucls:
#                         p = int(nucl == state)
#                         site_data.append(p)
#                     data.append(site_data)
#             aln_lens[part] = len(seq)

#         if ngenes > 1:
#             # get max set of parts. Only for multiple files
#             for node, parts in history.items():
#                 if len(parts) == ngenes:
#                     full_parts = parts.copy()
#                     break

#             # fill missing genes by '-'. Only for multiple files
#             for node, parts in history.items():
#                 if len(parts) != ngenes:
#                     unseen_parts = set(full_parts).difference(parts)
#                     for unp in unseen_parts:
#                         print(f"Gap filling for node {node}, part {unp}...", file=sys.stderr)
#                         for site in range(1, aln_lens[unp] + 1):
#                             site_data = [node, unp, site, "-", 0, 0, 0, 0]
#                             data.append(site_data)

#         df = pd.DataFrame(data, columns=columns)
#         return df

#     @staticmethod
#     def read_rates(path: str):
#         """ TODO write for many genes """
#         category = read_rates(path)
#         category = {"1": category}  # for each aln part
#         return category
    
#     @staticmethod
#     def get_mask(category: dict, cat_cutoff=1):
#         mask = {part: np.array(cats) >= cat_cutoff for part, cats in category.items()}
#         return mask

#     def close_db(self):
#         self.con.close()


def read_rates(path: str):
    """ TODO write for many genes """
    df = pd.read_csv(path, sep="\t", comment="#").sort_values("Site")
    category = df.Cat.values
    return category


def read_alignment(files: list, fmt="fasta"):
    """
    Read files alignments and prepare states table

    Arguments
    ---------
    files: list or iterable of string
        files with alignments; if some records not presented in all files table will be filled by '-' (gaps) 
    fmt: string
        format of files aligment; supported any format from biopython (fasta, phylip, etc.)
    
    Return
    ---------
        states: pd.DataFrame
    """
    ngenes = len(files)  
    columns = "Node Part Site State p_A p_C p_G p_T".split()
    nucls = "ACGT"
    aln_lens = dict()
    files = set(files)
    history = defaultdict(list)
    data = []
    for filepath in files:
        part = "1" if ngenes == 1 else os.path.basename(filepath).replace(".fna", "")
        alignment = SeqIO.parse(filepath, fmt)
        for rec in alignment:
            node = rec.name
            history[node].append(part)
            seq = str(rec.seq)
            for site, state in enumerate(seq, 1):
                site_data = [node, part, site, state]
                for nucl in nucls:
                    p = int(nucl == state)
                    site_data.append(p)
                data.append(site_data)
        aln_lens[part] = len(seq)

    if ngenes > 1:
        # get max set of parts. Only for multiple files
        for node, parts in history.items():
            if len(parts) == ngenes:
                full_parts = parts.copy()
                break

        # fill missing genes by '-'. Only for multiple files
        for node, parts in history.items():
            if len(parts) != ngenes:
                unseen_parts = set(full_parts).difference(parts)
                for unp in unseen_parts:
                    print(f"Gap filling for node {node}, part {unp}...", file=sys.stderr)
                    for site in range(1, aln_lens[unp] + 1):
                        site_data = [node, unp, site, "-", 0, 0, 0, 0]
                        data.append(site_data)

    df = pd.DataFrame(data, columns=columns)
    return df