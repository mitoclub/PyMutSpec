# CHANGELOG

<!-- ## 0.0.XX (2025-XX-XX)

Features:

- 

Fixes:

-  -->

## 0.0.15 (2026-07-16)

Features:

- Replaced `ete3` runtime dependency with a custom `TreeNode`/`Tree` implementation in `src/pymutspec/annotation/phylo_tree.py`
- Added BioPython-based Newick parsing and an ete3-compatible tree API (node iteration, branch traversal, distance computation, node search, and Newick serialisation)
- Added `tests/test_tree_vs_ete3.py` to verify parity with ete3 `PhyloTree` for tree loading and node/edge iteration

Fixes:

- Kept `ete3` only in `dev` extras for comparison tests
- Fixed `get_tree_len` name to `get_tree_height` to better reflect its purpose (computing tree height from root to leaves)
- Removed requirement for trees have outgroup and be named "ROOT"; now any tree can be used, and the ingroup is defined as the sister node of the outgroup if present, or the root if not.

**Full Changelog**: https://github.com/mitoclub/PyMutSpec/compare/0.0.14...0.0.15

## 0.0.14 (2025-03-31)

Features:

- Added new functions for parallel mutations extraction
- Outgroup now is optional node in the tree; ingroup now is just the sister node of outgroup and outgroup is the node that grow up from the root

Fixes:

- Repair package; pyproject.toml is work now. Setup.py deprecated, but simplified version saved.
- circular imports fixed
- all scripts now have main func

## 0.0.11 (2024-10-24)

Features:

- Edges sampling by @kpotoh in https://github.com/mitoclub/PyMutSpec/pull/9
- added barplots with errorbars
- citation added to readme

Fixes:

- fixed bug that allow short sequences be used in expected freqs calculation

**Full Changelog**: https://github.com/mitoclub/PyMutSpec/compare/0.0.10...0.0.11

## 0.0.10 (2024-04-14)

Features:

- added availability to collect non-syn mutational spectrum
- added a couple of arguments to calculate_mutspec func: scale, drop_underrepresented (this will drop mut types less than 0.9 by default)

Fixes:

- drop Times New Roman from default font for spectrum barplots
- deleted legacy log/doc files from repo
- added few tests for expected mutations collecting funcs
- fix scripts/calculate_mutspec.py : previously it failed when there was no mutations on branches, and reduced mutnum thresholds

## 0.0.8 (2023-11-13)

Features:

- Added little tests for plot functions
- Added few new functions from analysis notebooks, that used to tree spectra analysis

Fixes:

- Fixed uncertainty coefficient (phylocoef) calculation: based on only ingroup and geometric mean of distances from root to leaves
- General test for collect_mutations.py main class
