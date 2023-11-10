# PyMutSpec

Python library for mutational spectra analysis

<!-- Calculate mutational spectra using ancestral states from phylogenetic tree -->

## Requirements

- python 3.8+

## Installation

```bash
pip3 install pymutspec
```

https://pypi.org/project/PyMutSpec/

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

### Spectra barplots

<img src="./figures/ms12syn.png" width="300"/>

<img src="./figures/ms192syn.png" width="600"/>

## Links

1. [Iqtree](http://www.iqtree.org/) - efficient software for phylogenomic inference
2. [Genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1)
