import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from mutspec_utils.draw import plot_mutspec192kk, plot_mutspec192
import mutspec_utils

df = pd.read_csv("./data/processed/nematoda/dif_approaches/simple/mutspec192.tsv", sep="\t")
df = df[df.Label == "syn"]
print(df.shape)
print(df.head())

# plot_mutspec192kk(df, )
plot_mutspec192(df)
