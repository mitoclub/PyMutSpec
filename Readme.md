# Calculate mutational spectra using ancestral states from phylogenetic tree

## Environment

- python 3.8+

Activate python venv

```bash
python3 -m venv env_birds
source env_birds/bin/activate
pip install -r requirements.txt
```

## Workflow

6.1 Prepare appropriate format of leaves states

```bash
code...
```

## Stuff

...

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
