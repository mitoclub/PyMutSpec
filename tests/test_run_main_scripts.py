import os
import tmp

from scripts.collect_mutations import MutSpec

indir  = './tests/data/input'
outdir = './tests/data/output'
tree   = os.path.join(indir, 'gtr_100_cytb_replica_1.nwk')
states = os.path.join(indir, 'gtr_100_cytb_replica_1.fa')
rates  = os.path.join(indir, 'gtr_100_cytb_replica_1.inv')


# def test_collect_mutations_simple_cmd():
#     cmd = f'collect_mutations.py --tree {tree} --states {states} --states-fmt fasta --outdir {outdir} --gencode 2 --syn -q'


def test_mutspec_cls_basic(tmpdir):
    ms = MutSpec(tree, [states], tmpdir, 2, states_fmt='fasta', use_proba=False, syn=True)
    assert ms.MUT_LABELS == ['all', 'syn']
    rnd_genome = list(ms.get_random_genome().values())[0]
    assert len(rnd_genome) > 1000
    assert not ms.use_phylocoef

    ms.extract_mutspec_from_tree()


# def test_mutspec_cls_phylocoef(tmpdir):
#     ms = MutSpec(tree, [states], tmpdir, 2, states_fmt='fasta', use_proba=False, syn=True)
#     assert ms.MUT_LABELS == ['all', 'syn']
#     rnd_genome = list(ms.get_random_genome().values())[0]
#     assert len(rnd_genome) > 1000
#     assert not ms.use_phylocoef

#     ms.extract_mutspec_from_tree()
