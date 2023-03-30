# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2023-03-30
### Added
- The output file Base_report.cafe or Gamma_report.cafe has been added. This file has a 
format roughly the same as the CAFE 4 output, so tools and scripts that work with this
format should run properly.
- The tutorial now includes a cafe5_draw_tree.py script, demonstrating how to use BioPython
to create a tree.

### Changed
- Ordering of nodes has changed. The order should now be the same, relative to the source Newick
tree, as in the R package Ape.
- The model_change.tab file no longer includes a "+" sign before positive change counts.
- The simulator handles different max family and max root sizes correctly.
- The DocTest unit tester has been updated to v2.4.9.
- The configure script now checks for the existence of mkl.h before attempting to use it.
- Building the application should no longer fail if OpenMP is missing. (It will probably run
much slower, however)

### Removed
- Nothing


