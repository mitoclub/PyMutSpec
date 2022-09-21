## The “per-gene” test

For a single input tree and a large enough k, we have a distribution over signature values. Since we have limited data in this scenario, we use a parametric approach and fit a log-normal distribution to the signatures. Given a false positive tolerance rate α, we define values with a CDF above 1−α as outliers. Then, species associated with the outlier signatures are removed.

## The “all-gene” test

When a dataset includes several gene trees, all related by a species tree, combining the distributions across genes can increase the power. With many genes, we may also be able to distinguish outgroup species from outliers. The signatures of outgroups across all gene trees should be consistently higher than those of ingroups, and these high signatures will appear as part of the combined signature distribution. Thus, we may be able to avoid designating outgroups signatures as outliers.

In this test, we put the signature of all genes together to create one distribution. Unlike the per-gene test, here we have many data points, which enables us to use a non-parametric approach. We compute a kernel density function [28] over the empirical distribution of the combined set of signature values. To estimate the density, we use Gaussian kernels with Silverman’s rule of thumb smoothing bandwidth [28] (as implemented in the R package [29]). Given the density function and a false positive tolerance rate α, we define values with a CDF above 1−α as outliers.

## The “per-species” test

Outgroups can contribute to the tree diameter as much as erroneous species (Fig. 1b). To better distinguish outgroups from errors, when a set of gene trees are available, we can learn a distribution per species. Given a sufficient number of gene trees, the signatures of a species across all genes form a distribution that specifically captures the impact of that species on the gene tree diameter. These species-specific distributions naturally model the inherent difference between outgroups and ingroups in terms of their impacts on the tree diameter. More broadly, changes in the evolutionary tempo are captured naturally by the per species distributions.

In this test, we first compute a non-parametric distribution of the signature values for each species. When the signature of a species is not defined for a gene, we simply use zero as its signature. Then, for each species, we use the same non-parametric approach as in the all-gene test to compute a threshold for the signature value corresponding to the chosen α. Finally, we remove each species from those genes where its signature is strictly above its species-specific threshold.

## The default parameters of TreeShrink

TreeShrink has two parameters: α and k. By default, we set α to 0.05 (but users can choose other thresholds). Large values of k do not fit our goal of finding outlier species and can even lead to misleading results (e.g., Figure S2 in Appendix B, Additional file 1), but a small value of k may also miss outliers and may lead to insufficient data points for learning distributions.

Using a value of k that grows sublinearly with n (i.e., the number of leaves) gives us an algorithm that is fast enough for large n. For example, using 𝑘=Θ(𝑛√)
gives O(nh) running time, which on average is close to 𝑂(𝑛log𝑛) and is O(n2) in the worst case. While the choice must be ultimately made by the user, as a default, we set 𝑘=𝑚𝑖𝑛(𝑛4,5𝑛√). This heuristic formula ensures that our running time does not grow worse than quadratically with n but also avoids setting k to values close to n (thus also limits the proportion of leaves that could be removed).