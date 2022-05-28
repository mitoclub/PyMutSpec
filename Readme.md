# Calculate mutational spectra using ancestral states from phylogenetic tree

## Environment

- python 3.8+

Activate python venv

```bash
python3 -m venv env_birds
source env_birds/bin/activate
pip install -r requirements.txt
```

### Installing newick-utils-1.6

```bash
wget http://bioinfodbs.kantiana.ru/newick-utils-1.6.tar.gz
tar -xvzf newick-utils-1.6.tar.gz
cd newick-utils-1.6
./configure --prefix=/opt/newick-utils-1.6
make install
cd build/bin
ls
```

## Workflow

6.0 Root tree

```bash
nw_reroot anc.treefile Node3 > anc.treefile.rooted
```

6.1 Prepare appropriate format of leaves states

```bash
python scripts/6.1.terminal_genomes2iqtree_format.py --aln data/interim/alignments_birds_clean_clean --scheme data/interim/scheme_birds_genes.nex --out data/interim/leaves_birds_states.tsv
```

6.2 Prepare appropriate format of internal states from iqtree

```bash
python scripts/6.2.states2iqtree_format.py --anc data/interim/iqtree_runs/brun3/anc_kg/anc_kg.state --leaves data/interim/leaves_birds_states.tsv --out data/interim/anc_kg_states_birds.tsv
```

6.3 Caclulate MutSpec on most probable states

```bash
python scripts/6.3.calculate_mutational_spectra.py
```

6.4 Caclulate MutSpec using state probabilities

```bash
python scripts/6.4.calculate_mutational_spectra_proba.py
```

## Birds workflow

```bash
python mutspec/1.terminal_genomes2iqtree_format.py --aln data/example_birds/aln --scheme data/example_birds/scheme_birds_genes.nex --out data/example_birds/leaves_birds_states.tsv
python mutspec/2.states2iqtree_format.py --anc data/example_birds/anc_kg.state --leaves data/example_birds/leaves_birds_states.tsv --out data/example_birds/genes_states.tsv
# simple mutspec without probabilities

```

## Nematoda workflow

```bash
python mutspec/1.terminal_genomes2iqtree_format.py --aln data/example_nematoda/alignments_nematoda_clean --scheme data/example_nematoda/scheme_devilworm.nex --out data/example_nematoda/leaves_states_nematoda.tsv
python mutspec/2.states2iqtree_format.py --anc data/example_nematoda/anc_kg.state --leaves data/example_nematoda/leaves_states_nematoda.tsv --out data/example_nematoda/genes_states.tsv
# simple mutspec without probabilities
python mutspec/3.calculate_mutspec.py --tree data/example_nematoda/anc.treefile --anc data/example_nematoda/genes_states.tsv --leaves data/example_nematoda/leaves_states_nematoda.tsv
```

## Stuff

Plot trees

```bash
nw_display -s -S -v 30 -n -5 -b 'opacity:0' -i 'font-size:10px' -l 'font-family:serif;font-style:italic;font-size:small' -w 1500 anc.treefile.rooted > tree.svg
```

Plot trees

```bash
cd data/share
for fp in *.map; do sbs=`basename $fp .css.map`; nw_display -s -S -c $sbs.css.map -v 20 -b 'opacity:0' -i 'visibility:hidden' -l 'font-family:serif;font-style:italic;font-size:large' -d 'stroke-width:2' -w 1600 -R 30 anc.treefile.rooted > trees/tree_${sbs}.svg; done
```

## PastML

Used for calc mutspec without using custom phylo score

model подобрать

## Bottleneck

mode='db'

```txt
------------------------------
Function: extract_mutspec_from_tree
------------------------------
         11763856 function calls (11728284 primitive calls) in 32.786 seconds

   Ordered by: cumulative time, internal time
   List reduced from 1047 to 20 due to restriction <20>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.026    0.026   32.787   32.787 4.calculate_mutspec_proba.py:100(extract_mutspec_from_tree)
   783141    8.968    0.000   15.543    0.000 4.calculate_mutspec_proba.py:292(sample_context_fast)
       48    0.942    0.020   12.804    0.267 4.calculate_mutspec_proba.py:208(extract_mutations)
        8    0.236    0.029   11.612    1.452 utils.py:309(get_genome)
        8   11.326    1.416   11.326    1.416 {method 'execute' of 'sqlite3.Cursor' objects}
       48    1.215    0.025    5.532    0.115 4.calculate_mutspec_proba.py:264(collect_state_freqs)
      170    0.023    0.000    1.864    0.011 utils.py:224(calculate_mutspec)
```

mode='dict'

```txt
------------------------------
Function: extract_mutspec_from_tree
------------------------------
         11672319 function calls (11636779 primitive calls) in 22.308 seconds

   Ordered by: cumulative time, internal time
   List reduced from 1001 to 20 due to restriction <20>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.027    0.027   22.308   22.308 /home/mr/mutspec/mutspec/4.calculate_mutspec_proba.py:100(extract_mutspec_from_tree)
   783141    9.501    0.000   16.400    0.000 /home/mr/mutspec/mutspec/4.calculate_mutspec_proba.py:292(sample_context_fast)
       48    0.954    0.020   13.402    0.279 /home/mr/mutspec/mutspec/4.calculate_mutspec_proba.py:208(extract_mutations)
       48    1.280    0.027    5.860    0.122 /home/mr/mutspec/mutspec/4.calculate_mutspec_proba.py:264(collect_state_freqs)
      170    0.024    0.000    2.014    0.012 /home/mr/mutspec/mutspec/utils.py:224(calculate_mutspec)
      174    0.001    0.000    0.513    0.003 /home/mr/mutspec/mutspec/4.calculate_mutspec_proba.py:324(dump_table)
      174    0.002    0.000    0.502    0.003 /home/mr/env_bio/lib/python3.8/site-packages/pandas/core/generic.py:3388(to_csv)
```

## References

1. [Iqtree](http://www.iqtree.org/) - efficient software for phylogenomic inference
2. [Genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1)
