#!/usr/bin/env python3

import os
import sys
import re
from functools import partial

import click
import pandas as pd

from pymutspec.annotation import calculate_mutspec, rev_comp, transcriptor
from pymutspec.draw import plot_mutspec12, plot_mutspec192
from pymutspec.constants import possible_sbs12_set, possible_sbs192_set

ext_available = ['pdf', 'png', 'jpg']

@click.command("MutSpec visualizer", help="Visualize mutational spectra")
@click.option("-s", "--spectra", type=click.Path(True),  help="Path to spectra table, must contain columns Mut, ...")
@click.option("-o", '--outfile', type=click.Path(writable=True), default=None, help="Path to output image file, available extensions: ['pdf', 'png', 'jpg']")
@click.option("-t", "--invert", is_flag=True, help="Invert mutation (reverse complement), i.e. G>A --> C>T and A[C>T]G --> C[G>A]T")
def main(spectra, outfile, invert):
    ms = pd.read_csv(spectra, sep="\t")
    if outfile:
        ext = re.search(".+(\.\w+)$", outfile).group(1)
        if ext[1:] not in ext_available:
            raise ValueError("Extension for output file ({}) are not available, ({})".format(ext, ext_available))
        path_to_figure = outfile.replace(ext, "_{}" + ext)
    else:
        path_to_figure = spectra.replace(".tsv", "_{}.pdf")

    sbs_uniq = ms.Mut.unique()
    if possible_sbs12_set == set(sbs_uniq):
        plot_mutspec_func = partial(plot_mutspec12, style="bar")
        if invert:
            ms["Mut"] = ms["Mut"].str.translate(transcriptor)
    elif possible_sbs192_set == set(sbs_uniq):
        plot_mutspec_func = plot_mutspec192
        ms["Mut"] = ms["Mut"].apply(rev_comp)
    else:
        raise ValueError("Mut column must contain rates for all sbs (12 or 192)")

    if "Label" not in ms.columns:
        ms["Label"] = "all"
    for lbl, gr in ms.groupby("Label"):
        plot_mutspec_func(
            gr, 
            title=f"{lbl} mutational spectrum", 
            savepath=path_to_figure.format(lbl), 
            show=False,
        )


if __name__ == "__main__":
    main()
    # main("-b tmp/observed_mutations_iqtree.tsv -e tmp/expected_mutations.tsv -o tmp/ -l debug".split())
