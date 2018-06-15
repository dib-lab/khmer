# Change Log
All notable changes to the khmer project will be documented in this file.
See [keepachangelog](http://keepachangelog.com/) for more info.

The khmer project's command line scripts adhere to
[Semantic Versioning](http://semver.org/). The Python and C++ APIs are not yet
under semantic versioning, but will be in future versions of khmer.

## [Unreleased]

### Added
- Cython wrapper for liboxli.
- Cython containers for all liboxli classes.
- Header install for liboxli.
- New storage class using a Counting Quotient Filter with improved cache
  locality over bloom filters.
- New variants of the sequence bulk loading method with a "banding" mode and a
  "mask" mode. In "banding" mode, only k-mers whose hashed values fall within a
  specified range are counted. In "mask" mode, only k-mers not already pressent
  in the specified mask are counted.
    - `consume_seqfile_banding`
    - `consume_seqfile_with_mask`
    - `consume_seqfile_banding_with_mask`

### Changed
- Non-ACTG handling significantly changed so that only bulk-loading functions
  "clean" sequences of non-DNA characters. See #1590 for details.
- Split CPython wrapper file into per-class files under `src/khmer` and
  `include/khmer`.
- Moved liboxli headers to include/oxli and implementations to src/oxli.
- Removed all CPython wrappers except ReadParser and the standalone functions.
- Dropped support for Python 2.
- Changed to absolute imports.
- Some methods on LabelHash and Hashgraph have been changed to properties,
  or generators where appropriate.
- All constructors have been removed from khmer/__init__.py.
- GraphLabels does not inherit from Hashgraph.
- `trim-low-abund.py` doesn't error out when given multiple files with identical basenames

## [2.1.1] - 2017-05-25
### Added
- Document for submission to the Journal of Open Source Software.

### Fixed
- Several typos and outdated content in the documentation.

## [2.1] - 2017-05-21
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
- Several simple examples of the Python API and the C++ API in
  `examples/python-api` and `examples/c++-api`, respectively.
- New `assemble_linear_path` function for baiting unambiguous contigs with a
  seed k-mer from a hashtable.
- Support for assembling directly from k-mer graphs, and a new
  JunctionCountAssembler class.
- Add --info flag for obtaining citation information.
- Added Counttable and Nodetable to support non-reversible hashing
  functionality and k > 32.
- Add a new storage class using half a byte per entry. Exposed as
  SmallCounttable and SmallCountgraph.
- Added `cleaned_seq` attribute to `khmer.Read` class which provides a cleaned
  version of the sequence of each read.
- Added --summary-info to trim-low-abund.py to record run information in a file.
- `Nodetable`, `Counttable` and `SmallCounttable` use murmur hash 3 as hash
  function. This means they support kmers longer than 32 bases but means
  the hashes are not reversible.

### Changed
- Suppress display of -x and -N command line options in script help messages.
- Switch from nose to py.test as the testing framework.
- Switch from internally managed Jenkins setup to Travis CI for continuous
  integration testing.
- Renamed core data structures: CountingHash --> Countgraph,
  Hashbits --> Nodegraph.
- Replaced the IParser and FastxParser classes with a single ReadParser class.
  Different input formats are supported by templating ReadParser with a reader
  class.
- Renamed `consume_fasta` and related functions to `consume_seqfile`, with
  support for reading sequences from additional formats pending.
- Changed Sphinx documentation theme to "guzzle".

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
- Bug in `load-graph.py` that calculated required graph space incorrectly.
- Fix loading of empty partion map files

## [2.0] - 2015-10-08

Previous to the khmer 2.1 release, all changes were documented in a file named
`ChangeLog`. This file is now at `legacy/ChangeLog` for posterity.
