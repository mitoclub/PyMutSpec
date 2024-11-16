# PyMutSpec

[![PyPI Latest Release](https://img.shields.io/pypi/v/pymutspec.svg)](https://pypi.org/project/pymutspec/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/pymutspec.svg?label=PyPI%20downloads)](
https://pypi.org/project/pymutspec/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pymutspec)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/mitoclub/pymutspec)
![GitHub last commit](https://img.shields.io/github/last-commit/mitoclub/pymutspec)

[![NeMu Paper](https://img.shields.io/badge/DOI-10.1093%2Fnar%2Fgkae438-blue)](https://doi.org/10.1093/nar/gkae438)

<!-- [![CI - Test](https://github.com/pandas-dev/pandas/actions/workflows/unit-tests.yml/badge.svg)](https://github.com/pandas-dev/pandas/actions/workflows/unit-tests.yml) -->
<!-- [![Coverage](https://codecov.io/github/pandas-dev/pandas/coverage.svg?branch=main)](https://codecov.io/gh/pandas-dev/pandas) -->

Python library for mutational spectra analysis

<!-- Calculate mutational spectra using ancestral states from phylogenetic tree -->

## Installation

```bash
pip3 install pymutspec
```

## Example code

```python
from Bio import SeqIO
from pymutspec.annotation import calculate_mutspec, CodonAnnotation
from pymutspec.draw import plot_mutspec12, plot_mutspec192

coda = CodonAnnotation(gencode=2) # mitochondrial genetic code

path_to_observed_mutations = ... 
path_to_reference_seq = ...

# load data (mutations and sequence)
gene = SeqIO.parse(path_to_reference_seq, format='fasta')
observed_mutations = pd.read_csv(path_to_observed_mutations, sep='\t')
for col in ['Mut', 'MutType']:
    assert col in observed_mutations.columns

# sample only syn mutations
mut_syn = observed_mutations[observed_mutations.MutType >= 1] # 0 for all mutations, 1 for syn, 2 for fourfold syn (syn4f)

# derive expected mutations from reference gene
sbs12_freqs, sbs192_freqs = coda.collect_exp_mut_freqs(gene, labels['all', 'syn', 'syn4f'])
sbs12_freqs_syn = sbs12_freqs['syn']
sbs192_freqs_syn = sbs192_freqs['syn']

# calculate mutation spectra
spectra12 = calculate_mutspec(mut_syn, sbs12_freqs_syn, use_context=False)
spectra192 = calculate_mutspec(mut_syn, sbs192_freqs_syn, use_context=True)

# plot mutation spectra
plot_mutspec12(spectra12)
plot_mutspec192(spectra192)
```

### Example spectra barplots

<img src="https://raw.githubusercontent.com/mitoclub/PyMutSpec/master/figures/ms12syn.png" width="300"/>

<img src="https://raw.githubusercontent.com/mitoclub/PyMutSpec/master/figures/ms192syn.png" width="600"/>

## Links

1. [IQ-Tree2](http://www.iqtree.org/) - efficient software for phylogenomic inference
2. [Genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1)

## How to cite?

    Bogdan Efimenko, Konstantin Popadin, Konstantin Gunbin, NeMu: a comprehensive pipeline for accurate reconstruction of neutral mutation spectra from evolutionary data, Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W108–W115, https://doi.org/10.1093/nar/gkae438
