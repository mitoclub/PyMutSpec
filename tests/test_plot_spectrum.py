import pandas as pd

from pymutspec.draw import plot_mutspec12, plot_mutspec192
from pymutspec.draw.sbs_orders import ordered_sbs192_kk


def test_plot_mutspec12():
    ms12 = pd.read_csv('./tests/data/ms12syn_iqtree.tsv', sep='\t')
    ax1 = plot_mutspec12(ms12, style='bar')
    ax2 = plot_mutspec12(ms12, style='box')


def test_plot_mutspec192():
    ms192 = pd.read_csv('./tests/data/ms192syn_iqtree.tsv', sep='\t')
    ax1 = plot_mutspec192(ms192, style='bar', labels_style="cosmic",)
    ax2 = plot_mutspec192(ms192, style='box', labels_style="long",)
    ax3 = plot_mutspec192(ms192, sbs_order=ordered_sbs192_kk)

    ms192 = ms192[(ms192.ObsNum > 1) & (ms192.ExpNum > 1)]
    ax4 = plot_mutspec192(ms192)
