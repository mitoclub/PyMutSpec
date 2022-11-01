#!/bin/bash

cd /home/kpotoh/mutspec-utils/tmp/evolve_sc

raw_tree=../evolve/iqtree_anc_tree.nwk
spectra=../evolve/ms12syn_iqtree.tsv
mulal=../evolve/alignment_checked.fasta

label=iqtree
replics=10
GENCODE=2
scale_tree=1

tree=tree.nwk
nw_prune $raw_tree OUTGRP > $tree
python3 ../../scripts/pyvolve_process.py -a $mulal -t $tree -s $spectra -o seqfile.fasta -r $replics --write_anc -c $GENCODE -l $scale_tree

echo > pyvolve_full_$label.log

for fasta_file in seqfile_sample-*.fasta; do
	echo $fasta_file
	python3 ../../scripts/alignment2iqtree_states.py $fasta_file $fasta_file.state
	python3 ../../scripts/3.collect_mutations.py --tree $tree --states $fasta_file.state --gencode $GENCODE --syn --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_full_$label.log
	echo -e "\n\n">> pyvolve_full_$label.log
	cat mout/mutations.tsv > $fasta_file.mutations
done

python3 ../../scripts/concat_mutations.py seqfile_sample-*.fasta.mutations full_mutations.txt

mkdir -p out

# TODO add additional agrs
python3 ../../scripts/calculate_mutspec_pyvolve.py -b full_mutations.txt -e mout/expected_mutations.tsv \
	-o out -l iqtree --outgrp OUTGRP


cd /home/kpotoh/mutspec-utils
