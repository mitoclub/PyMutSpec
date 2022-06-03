import os
import re
import sys
import sqlite3
from collections import defaultdict
from typing import Dict, List

import numpy as np
import pandas as pd
from Bio import SeqIO
import tqdm


def load_scheme(path: str) -> Dict[str, str]:
    """
    parse files like scheme_birds_genes.nex (just separated genes)

    return dict(charset_lbl: gene_fp)
    """
    with open(path) as handle:
        raw_file = handle.read()
    charsets = re.findall("charset\s(\w+)\s?=\s?([\w_\.]+)(\s?:.+)?;", raw_file)
    scheme = {i: os.path.basename(fp) for i, (_, fp, _) in enumerate(charsets, 1)}
    return scheme


def get_aln_files(path: str):
    assert os.path.isdir(path), "path is not directory"
    raw_files = os.listdir(path)
    files = set(
        [os.path.join(path, x) for x in raw_files if x.endswith(".fna")]
    )
    return files
    

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


def read_genbank_ref(path: str):
    gb_file = next(SeqIO.parse(path, "genbank"))
    ftypes_nc = {'rRNA', 'tRNA'}
    full_nucls = set("ACGT")
    data = []
    df: pd.DataFrame = None
    for ftr in gb_file.features:
        if ftr.type == "source":
            source = ftr.extract(gb_file)
            seq = str(source.seq)
            for pos, nuc in enumerate(seq):
                context = seq[pos - 1: pos + 2]
                if len(context) < 3 or len(set(context).difference(full_nucls)) != 0:
                    context = None
                if nuc not in full_nucls:
                    nuc = context = None
                data.append({"Pos": pos + 1, "Nuc": nuc, "Context": context})
            df = pd.DataFrame(data)
            df["Strand"] = 0
            continue

        for pos in list(ftr.location):
            df.at[pos, "Type"] = ftr.type
            df.at[pos, "Strand"] = ftr.strand        
            if ftr.type == 'CDS' or ftr.type in ftypes_nc:
                df.at[pos, "GeneName"] = ftr.qualifiers["gene"][0]

    # add codon features
    df["PosInGene"] = -1
    df["PosInCodon"] = -1
    for gene_name in df[(df.Type == "CDS") & (df.Strand == 1)].GeneName.unique():
        gdf = df[df.GeneName == gene_name]
        seq = gdf.Nuc.values
        for pos_in_gene, pos in enumerate(gdf.index):
            pic = pos_in_gene % 3
            codon = seq[pos_in_gene - pic: pos_in_gene - pic + 3]
            codon = "".join(codon) if len(set(codon).difference(full_nucls)) == 0 else None
            df.at[pos, "Codon"] = codon
            df.at[pos, "PosInGene"] = pos_in_gene + 1
            df.at[pos, "PosInCodon"] = pic + 1

    df["Strand"] = df["Strand"].astype(np.int8)
    df["PosInCodon"] = df["PosInCodon"].astype(np.int8)
    df["PosInGene"] = df["PosInGene"].astype(np.int32)
    return df
