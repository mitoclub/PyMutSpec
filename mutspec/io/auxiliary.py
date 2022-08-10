import os
import re
from typing import Dict

import numpy as np
import pandas as pd
from Bio import SeqIO


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
            cdn = seq[pos_in_gene - pic: pos_in_gene - pic + 3]
            cdn = "".join(cdn) if len(set(cdn).difference(full_nucls)) == 0 else None
            df.at[pos, "Codon"] = cdn
            df.at[pos, "PosInGene"] = pos_in_gene + 1
            df.at[pos, "PosInCodon"] = pic + 1

    df["Strand"] = df["Strand"].astype(np.int8)
    df["PosInCodon"] = df["PosInCodon"].astype(np.int8)
    df["PosInGene"] = df["PosInGene"].astype(np.int32)
    return df
