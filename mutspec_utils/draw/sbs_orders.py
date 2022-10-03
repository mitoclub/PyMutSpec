import pandas as pd

from mutspec_utils.constants import possible_sbs192
from mutspec_utils.annotation import rev_comp, translator

kk_lbls = "A>C A>G A>T C>T G>C G>T".split()
cosmic_lbls = "C>A C>G C>T T>A T>C T>G".split()

df = pd.DataFrame({"sbs": possible_sbs192})
df["sbs_base"] = df["sbs"].str.slice(2, 5)
df["sbs_base_revcomp"] = df["sbs_base"].str.translate(translator)
df["sbs_revcomp"] = df["sbs"].apply(rev_comp)
df["is_cosmic"] = df["sbs_base"].isin(cosmic_lbls)
df["is_kk"] = df["sbs_base"].isin(kk_lbls)
df["sbs_base_for_sorting_kp"] = df.apply(
    lambda x: x.sbs_base + "1" if x.is_cosmic else x.sbs_base_revcomp + "2", axis=1)
df["sbs_for_ordering_kk"] = df.apply(lambda x: x.sbs if x.is_kk else x.sbs_revcomp, axis=1)
df["sbs_for_ordering_kp"] = df.apply(lambda x: x.sbs if x.is_cosmic else x.sbs_revcomp, axis=1)

ordered_sbs192_kp = list(df.sort_values(["sbs_base_for_sorting_kp", "sbs_for_ordering_kp"]).sbs.values)
ordered_sbs192_kk = list(df.sort_values(["sbs_base", "sbs_for_ordering_kk"]).sbs.values)
del df
