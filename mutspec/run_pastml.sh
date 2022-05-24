# example
pastml -t data/example_pastml/Albanian.tree.152tax.tre -d data/example_pastml/data.txt --html data/example_pastml/tree.html --html_compressed data/example_pastml/map.html --data_sep ,

# birds
pastml -t data/example_birds/anc_kg.treefile -d data/leaves_birds.pastml.tsv --html_compressed data/example_birds/map.html --html data/example_birds/tree.html -v --threads 24 -m HKY # -m HKY

# birds
pastml -t data/example_nematoda/anc.treefile -d data/example_nematoda/leaves.pastml.tsv --html_compressed data/example_nematoda/map.html --html data/example_nematoda/tree.html -v --threads 24 # -m HKY


# -o gives combined_ancestral_states.tab
