# Change Log
All notable changes to the khmer project will be documented in this file.
See [keepachangelog](http://keepachangelog.com/) for more info.

The khmer project's command line scripts adhere to
[Semantic Versioning](http://semver.org/). The Python and C++ APIs are not yet
under semantic versioning, but will be in future versions of khmer.

## [Unreleased]
### Added
- New `--no-reformat` option for `interleave-reads.py` script disables default
  read name correction behavior.
- New `HashSet` data structure for managing collections of k-mer hashes and
  tags.
- khmer package version now included in `.info` files.
- New `-o|--outfile` option for `filter-abund-single.py` script.
- New sandbox script `extract-compact-dbg.py` for computing a compact de Bruijn
  graph from sequence data.
- New `--quiet` flag to several scripts, silencing diagnostic messages in
  terminal output.
- Support for human-friendly memory requests (2G instead of 2000000000 or 2e9).
- Support for variable-coverage trimming in the `filter-abund-single.py` script.
- A simple example of the C++ API in `examples/c++-api`.
- New `assemble_linear_path` function for baiting unambiguous contigs with a
  seed k-mer from a hashtable.
- Support for assembling directly from k-mer graphs, and a new
  JunctionCountAssembler class.
- Add --info flag for obtaining citation information.
- Add a new storage class using half a byte per entry. Exposed as
  SmallCounttable and SmallCountgraph.

### Changed
- Switch from nose to py.test as the testing framework.
- Switch from internally managed Jenkins setup to Travis CI for continuous
  integration testing.
- Renamed core data structures: CountingHash --> Countgraph,
  Hashbits --> Nodegraph.

### Fixed
- Bug in compressed(gzip) streaming output from scripts
- The hashbits `update_from` function to correctly track occupied bins for
  calculating FPR.
- Bug in the `filter-abund.py` script when `--gzip` and `-o` flags are used
  simultaneously.
- Bug in the hashtable `get_kmers` function based on incorrect usage of the
  `substr` function.
- Bug in `broken_paired_reader` related to dropping short reads when
  `require_paired` is set.
- Bug related to handling lowercase [acgtn] characters in input data.

## [2.0] - 2015-10-08

Previous to the khmer 2.1 release, all changes were documented in a file named
`ChangeLog`. This file is now at `legacy/ChangeLog` for posterity.
