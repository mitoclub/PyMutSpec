# PyMutSpec

Python library for mutational spectra analysis

<!-- Calculate mutational spectra using ancestral states from phylogenetic tree -->

## Install

- python 3.8+

Install python venv (optionally)

```bash
python3 -m venv env_spectra
source env_spectra/bin/activate
# pip install -r requirements.txt
python3 setup.py install
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

## References

1. [Iqtree](http://www.iqtree.org/) - efficient software for phylogenomic inference
2. [Genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1)
