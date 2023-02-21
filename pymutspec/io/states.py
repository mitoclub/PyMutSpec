import os
import sys
import sqlite3
from collections import defaultdict
from typing import List
from unicodedata import category

import numpy as np
import pandas as pd
import tqdm


class GenesStates:
    """
    tips:
    - use mode="dict" if your mtDNA tree is small (~1500 nodes x 10,000 nucleotides require ~350MB of RAM); 
    big trees with 100,000 nodes require tens of GB of RAM, therefore use mode="db"
    """
    def __init__(self, path_states: List[str], path_to_db=None, mode="dict", rewrite=False, use_proba=True, path_to_rates=None, cat_cutoff=2):
        self.mode = mode
        self.use_proba = use_proba
        print(f"Genome states storage mode = '{mode}'", file=sys.stderr)
        self.path_to_db = path_to_db
        if mode == "dict":
            self._prepare_node2genome(path_states)
        elif mode == "db":
            if path_to_db is None:
                raise ValueError("Pass the path to database to use or write it")
            else:
                self._prepare_db(path_states, rewrite)
        else:
            raise ValueError("Mode must be 'dict' or 'db'")
        
        genes_sizes = {part: len(x) for part, x in self.get_random_genome().items()}
        self.genome_size = sum(genes_sizes.values())
        self.category = self.read_rates(path_to_rates) if path_to_rates else None
        self.mask = None

        if self.category: 
            for part in genes_sizes:
                if genes_sizes[part] != len(self.category[part]):
                    raise RuntimeError("Wrong rates: number of positions in alignment are not equal to number of rates", file=sys.stderr)
            self.mask = self.get_mask(self.category, cat_cutoff)
    
    def get_genome(self, node: str):
        if self.mode == "dict":
            return self.node2genome[node]
            
        elif self.mode == "db":
            genome_raw = defaultdict(list)
            cur = self.con.cursor()
            if self.use_proba:
                query = f"""SELECT Part, Site, p_A, p_C, p_G, p_T FROM states WHERE Node='{node}'"""
                dtype = [
                    ("Site", np.int32),
                    ("p_A", np.float32), ("p_C", np.float32), 
                    ("p_G", np.float32), ("p_T", np.float32),
                ]
            else:
                query = f"""SELECT Part, Site, State FROM states WHERE Node='{node}'"""
                dtype = [("Site", np.int32), ("State", np.object_)]

            for row in cur.execute(query):
                part = row[0]
                state = row[1:]
                genome_raw[part].append(state)
            
            genome = {}
            for part, state_tup in genome_raw.items():
                state = np.array(state_tup, dtype=dtype)
                state_sorted = np.sort(state, order="Site")
                state_sorted_devoid = np.array(list(map(list, state_sorted)))
                state = state_sorted_devoid[:, 1:] if self.use_proba else state_sorted_devoid[:, 1]
                genome[part] = state
            return genome

        else:
            raise ValueError("mode must be 'dict' or 'db'")
    
    def get_random_genome(self):
        node = next(iter(self.nodes))
        return self.get_genome(node)

    def _prepare_db(self, path_states, rewrite=False):
        """sequentially read tsv and write to db"""
        done = False
        if os.path.exists(self.path_to_db):
            if rewrite:
                os.remove(self.path_to_db)
            else:
                print(f"Loading existing database from {self.path_to_db}", file=sys.stderr)
                done = True
        try:
            con = sqlite3.connect(self.path_to_db)
        except:
            raise ValueError("Cannot connect to database by given path")
        
        cur = con.cursor()
        self.con = con
        if done:
            self.nodes = self._get_nodes_from_db()
            return

        print("Writing new database file", file=sys.stderr)
        cur.execute('''CREATE TABLE states
               (Node TEXT, Part INTEGER, Site INTEGER, State TEXT, p_A REAL, p_C REAL, p_G REAL, p_T REAL)'''
        )
        for path in path_states:
            if path is None:
                continue
            print(f"Loading {path} ...", file=sys.stderr)
            handle = open(path, "r")
            header = handle.readline().strip()
            if header != "Node\tPart\tSite\tState\tp_A\tp_C\tp_G\tp_T":
                handle.close()
                con.close()
                raise ValueError(f"Inappropriate type of table, expected another columns order,\ngot {repr(header)}")

            for line in tqdm.tqdm(handle, total=8652300):  # TODO estimate total
                row = line.strip().split()
                query = "INSERT INTO states VALUES ('{}',{},{},'{}',{},{},{},{})".format(*row)
                cur.execute(query)

            con.commit()
            handle.close()

        self.nodes = self._get_nodes_from_db()
        # con.close()
    
    def _get_nodes_from_db(self):
        nodes = set()
        cur = self.con.cursor()
        for node in cur.execute("SELECT DISTINCT Node from states"):
            nodes.add(node[0])
        return nodes

    def _prepare_node2genome(self, path_states, states_dtype=np.float32):
        dtypes = {
            "p_A":  states_dtype, "p_C": states_dtype, 
            "p_G":  states_dtype, "p_T": states_dtype,
            "Site": np.int32, "Node": str,     #"Part": np.int8,
        }
        if self.use_proba:
            usecols = ["Node", "Part", "Site", "p_A", "p_C", "p_G", "p_T"]
        else:
            usecols = ["Node", "Part", "Site", "State"]
        node2genome = defaultdict(dict)
        for path in path_states:
            if path is None:
                continue
            print(f"Loading {path}...", file=sys.stderr)
            states = pd.read_csv(path, sep="\t", comment='#', usecols=usecols, dtype=dtypes)
            aln_sizes = states.groupby("Node").apply(len)
            assert aln_sizes.nunique() == 1, "uncomplete leaves state table"
            states = states.sort_values(["Node", "Part", "Site"])
            gr = states.groupby(["Node", "Part"])
            for (node, part), gene_pos_ids in tqdm.tqdm(gr.groups.items()):
                gene_df = states.loc[gene_pos_ids]  # TODO simplify: iterate over gr
                if self.use_proba:
                    gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
                else:
                    gene_states = gene_df.State.values
                node2genome[node][part] = gene_states

        self.node2genome = node2genome
        self.nodes = set(node2genome.keys())

    @staticmethod
    def read_rates(path: str):
        """ TODO write for many genes """
        category = read_rates(path)
        category = {1: category}  # for each aln part
        return category
    
    @staticmethod
    def get_mask(category: dict, cat_cutoff=2):
        mask = {part: np.array(cats) >= cat_cutoff for part, cats in category.items()}
        return mask

    def close_db(self):
        self.con.close()


def read_rates(path: str):
    """ TODO write for many genes """
    df = pd.read_csv(path, sep="\t", comment="#").sort_values("Site")
    category = df.Cat.values
    return category
