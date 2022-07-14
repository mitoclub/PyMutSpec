parallel iqtree2 -s {} -m MFP -nt 2 --prefix data/example_nematoda/trees_with_gaps/{/.} ::: data/example_nematoda/alignments_nematoda_clean/*
parallel iqtree2 -s {} -m MFP -nt 2 --prefix data/example_nematoda/trees/{/.} ::: data/example_nematoda/trimed_aln_nematoda/*
