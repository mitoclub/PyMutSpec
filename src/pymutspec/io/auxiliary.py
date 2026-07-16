import os
import re
from typing import Dict


def load_scheme(path: str) -> Dict[str, str]:
    """
    TODO deprecate
    parse files like scheme_birds_genes.nex (just separated genes)

    return dict(charset_lbl: gene_fp)
    """
    with open(path) as handle:
        raw_file = handle.read()
    charsets = re.findall(r"charset\s(\w+)\s?=\s?([\w_\.]+)(\s?:.+)?;", raw_file)
    scheme = {i: os.path.basename(fp) for i, (_, fp, _) in enumerate(charsets, 1)}
    return scheme


def get_aln_files(path: str):
    """
    Return the set of ``*.fna`` alignment files in a directory.

    Arguments
    ---------
    path: str
        Path to a directory containing alignment files.

    Return
    ------
    files: set of str
        Absolute paths to all ``*.fna`` files found in *path*.

    Raises
    ------
    AssertionError
        If *path* is not a directory.
    """
    assert os.path.isdir(path), "path is not directory"
    raw_files = os.listdir(path)
    files = set(
        [os.path.join(path, x) for x in raw_files if x.endswith(".fna")]
    )
    return files
