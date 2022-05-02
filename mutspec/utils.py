import cProfile
import pstats
import re
import sys
from collections import defaultdict
from typing import List, Set, Tuple, Union

from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA

PATH_TO_GENCODE5 = "data/external/genetic_code5.txt"


def _read_gencode(path: str) -> List[Tuple[str]]:
    pattern = re.compile("([ACGT]{3})\s([A-Z\*])\s([A-Za-z]{3})\s(i?)")
    gencode = []
    with open(path) as fin:
        for line in fin:
            for code in pattern.findall(line):
                gencode.append(code)
    return gencode


def _read_start_stop_codons(path: str) -> Tuple[Set[str]]:
    gencode = _read_gencode(path)
    startcodons = set()
    stopcodons = set()
    for code in gencode:
        if code[-1] == "i":
            startcodons.add(code[0])
        if code[1] == "*":
            stopcodons.add(code[0])
    return startcodons, stopcodons


def __prepare_codontable(codontable: Union[NCBICodonTableDNA, int]):
    if isinstance(codontable, NCBICodonTableDNA):
        pass
    elif isinstance(codontable, int):
        codontable = CodonTable.unambiguous_dna_by_id[codontable]
    else:
        ValueError("passed codontable is not appropriate")
    return codontable


def read_start_stop_codons(codontable: Union[NCBICodonTableDNA, int]):
    codontable = __prepare_codontable(codontable)
    return set(codontable.start_codons), set(codontable.stop_codons)


def _extract_ff_codons(codontable: Union[NCBICodonTableDNA, int]):
    codontable = __prepare_codontable(codontable)
    aa2codons = defaultdict(set)
    for codon, aa in codontable.forward_table.items():
        aa2codons[aa].add(codon)

    ff_codons = set()
    for aa, codons in aa2codons.items():
        if len(codons) >= 4:
            interim_dct = defaultdict(set)
            for codon in codons:
                interim_dct[codon[:2]].add(codon)

            for nn in interim_dct:
                if len(interim_dct[nn]) == 4:
                    ff_codons = ff_codons.union(interim_dct[nn])
    return ff_codons


def extract_syn_codons(codontable: Union[NCBICodonTableDNA, int]):
    """ extract synonymous (codons that mutate without amino acid change) 
    and fourfold codons from codon table

    usefull function for expected mutspec (filtration)

    return mapping[(cdn, pic)] of syn codons and set of ff codons
    """
    codontable = __prepare_codontable(codontable)
    aa2codons = defaultdict(set)
    for codon, aa in codontable.forward_table.items():
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


def is_syn_codons(codon1: str, codon2: str, codontable: Union[NCBICodonTableDNA, int]):
    """
    extract codons containing mutation that are synonymous

    return dict[codon: set[PosInCodon]]
    """
    codontable = __prepare_codontable(codontable)
    gc = codontable.forward_table
    return gc.get(codon1, "*") == gc.get(codon2, "*")


def node_parent(node):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None


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
    # print(read_start_stop_codons(2))
    # print(extract_syn_codons(2))
    code = 5
    syn_codons, ff_codons = extract_syn_codons(code)
    ff_codons_true = extract_ff_codons(code)
    print(syn_codons)
    assert ff_codons == ff_codons_true

