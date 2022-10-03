import os
import sys
import sqlite3
from collections import defaultdict
from typing import List

import numpy as np
import pandas as pd
from Bio import SeqIO
import tqdm


class GenomeStates:
    """
    tips:
    - use mode="dict" if your mtDNA tree is small (~1500 nodes require ~350MB of RAM); 
    big trees with 100000 nodes require tens GB of RAM, therefore use mode="db"
    """
    def __init__(self, path_states: List[str], path_to_db=None, mode="dict", rewrite=False, proba_mode=True):
        self.mode = mode
        self.proba_mode = proba_mode
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
        
        self.genome_size = sum([len(part) for _, part in self.get_random_genome().items()])
    
    def get_genome(self, node: str):
        if self.mode == "dict":
            return self.node2genome[node]
            
        elif self.mode == "db":
            genome_raw = defaultdict(list)
            cur = self.con.cursor()
            if self.proba_mode:
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
                state = state_sorted_devoid[:, 1:] if self.proba_mode else state_sorted_devoid[:, 1]
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
            self.nodes = self._get_nodes()
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

        self.nodes = self._get_nodes()
        # con.close()
    
    def _get_nodes(self):
        nodes = set()
        cur = self.con.cursor()
        for node in cur.execute("SELECT DISTINCT Node from states"):
            nodes.add(node[0])
        return nodes

    def _prepare_node2genome(self, path_states, states_dtype=np.float32):
        dtypes = {
            "p_A":  states_dtype, "p_C": states_dtype, 
            "p_G":  states_dtype, "p_T": states_dtype,
            "Site": np.int32,     #"Part": np.int8,
        }
        if self.proba_mode:
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
            for (node, part), gene_pos_ids in gr.groups.items():
                gene_df = states.loc[gene_pos_ids]
                if self.proba_mode:
                    gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
                else:
                    gene_states = gene_df.State.values
                node2genome[node][part] = gene_states

        self.node2genome = node2genome
        self.nodes = set(node2genome.keys())
    
    def close_db(self):
        self.con.close()
