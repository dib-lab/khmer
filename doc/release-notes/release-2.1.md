# khmer v2.1 release notes

We are pleased to announce release version 2.1 of khmer!
This release includes several new features, bug fixes, and internal changes.
[CHANGELOG.md](https://github.com/dib-lab/khmer/blob/v2.1/CHANGELOG.md) includes a complete description of changes made since the previous release.
A concise summary of the most relevant changes is provided below.

The latest version of the khmer documentation is available at https://khmer.readthedocs.org/en/v2.1/.


## Items of note

### New features

- New storage class using half a byte per entry. Exposed as SmallCounttable and SmallCountgraph.
- New Counttable, SmallCounttable, and Nodetable classes to support non-reversible hashing functionality and k > 32.
- Support for human-friendly memory requests (2G) in addition to the previous style of requests (2000000000 or 2e9).
- Support for variable-coverage trimming in the `filter-abund-single.py` script.
- khmer package version now included in `.info` files.
- Several simple examples of the Python API and the C++ API in `examples/python-api` and `examples/c++-api`, respectively.
- Support for assembling directly from k-mer graphs, and a new JunctionCountAssembler class.
- New sandbox script `extract-compact-dbg.py` for computing a compact de Bruijn graph from sequence data.

### Bug fixes

- Streaming gzip-compressed output from scripts now works correctly.
- The `load-graph.py` script now calculates required graph space correctly.
- The `broken_paired_reader` no longer drops short reads when `require_paired` is set.

### Other changes

- Command-line options `-x` and `-N` now suppressed by default in script help messages.
- Renamed core data structures: CountingHash --> Countgraph, Hashbits --> Nodegraph.
- Replaced the `IParser` and `FastxParser` classes with a single `ReadParser` class. Different input formats are supported by templating `ReadParser` with a reader class.
- Renamed `consume_fasta` and related functions to `consume_seqfile`, with support for reading sequences from additional formats pending.


## Known issues:

- `load-graph.py` in multithreaded mode will find slightly different number of unique kmers. This is being investigated in https://github.com/dib-lab/khmer/issues/1248
- Any scripts that consume FASTA or FASTQ data are unresponsive to attempts to terminate the program with `ctrl-c`. Eventually khmer should handle this properly, but for now it should be possible to halt a script using `ctrl-\`. This is being tracked in https://github.com/dib-lab/khmer/issues/1618.

## Contributors

‡@aaliyari, @ctb, ‡@caitsydney, @camillescott, @standage, @kdmurray91, @ljcohen, @luizirber, @mr-c, ‡@NBKingsley, ‡@ryneches, ‡@rcs333, ‡@satta, ‡@shannonekj, ‡@betatim

‡ Indicates new contributors
