#!/bin/bash

raw_tree=tmp/evolve/iqtree_anc_tree.nwk
spectra=tmp/ms12syn_.tsv
mulal=tmp/evolve/alignment_checked.fasta

label="iqtree"
replics=50
GENCODE=2
scale_tree=1

tree=tmp/evolve_sc/tree.nwk
nw_prune $raw_tree OUTGRP > $tree
python3 scripts/pyvolve_process.py -a $mulal -t $tree -s $spectra -o tmp/evolve_sc/seqfile.fasta -r $replics --write_anc -c $GENCODE -l $scale_tree

echo > tmp/evolve_sc/pyvolve_full_$label.log

for fasta_file in tmp/evolve_sc/seqfile_sample-*.fasta; do
	echo $fasta_file
	python3 scripts/alignment2iqtree_states.py $fasta_file $fasta_file.state
	python3 scripts/3.collect_mutations.py --tree $tree --states $fasta_file.state --gencode $GENCODE --syn --no-mutspec --outdir tmp/evolve_sc/mout --force
	cat tmp/evolve_sc/mout/run.log >> tmp/evolve_sc/pyvolve_full_$label.log
	cat tmp/evolve_sc/mout/mutations.tsv > $fasta_file.mutations
	
	echo >> tmp/evolve_sc/pyvolve_full_$label.log
	echo >> tmp/evolve_sc/pyvolve_full_$label.log
done

python3 scripts/concat_mutations.py tmp/evolve_sc/seqfile_sample-*.fasta.mutations tmp/evolve_sc/full_mutations.txt

mkdir -p tmp/evolve_sc/out
python3 scripts/calculate_mutspec_pyvolve.py -b tmp/evolve_sc/full_mutations.txt -e tmp/evolve_sc/mout/expected_mutations.tsv \
	-o tmp/evolve_sc/out -l iqtree --outgrp OUTGRP

# python3 scripts/calculate_mutspec_pyvolve.py -b tmp/evolve_sc/mout/mutations.tsv -e tmp/evolve_sc/mout/expected_mutations.tsv \
# 	-o . --outgrp OUTGRP -m $mnum192 -l $label