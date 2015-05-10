Sandbox scripts
===============

Scripts in this directory, 'sandbox', are various utility or trial
scripts that we have not fully tested.  They are also not under
semantic versioning, so their functionality and command line arguments
may change without notice.

We are still triaging and documenting the various scripts.

----

Awaiting promotion to scripts:

* calc-error-profile.py - calculate a per-base "error profile" for shotgun sequencing data, w/o a reference. (Used/tested in `2014 paper on semi-streaming algorithms <https://github.com/ged-lab/2014-streaming/blob/master/>`__)
* count-kmers.py - output k-mer counts for multiple input files.
* count-kmers-single.py - output k-mer counts for a single k-mer file.
* correct-errors.py - streaming error correction.
* unique-kmers.py - estimate the number of k-mers present in a file with the HyperLogLog low-memory probabilistic cardinality estimation algorithm.

Scripts with recipes:

* calc-median-distribution.py - plot coverage distribution; see `khmer-recipes #1 <https://github.com/dib-lab/khmer-recipes/tree/master/001-extract-reads-by-coverage>`__
* collect-reads.py - subsample reads until a particular average coverage; see `khmer-recipes #2 <https://github.com/dib-lab/khmer-recipes/tree/master/002-collect-subset-of-high-coverage>`__
* saturate-by-median.py - calculate collector's curve on shotgun sequencing; see `khmer-recipes #4 <https://github.com/dib-lab/khmer-recipes/tree/master/004-estimate-sequencing-saturation>`__
* slice-reads-by-coverage.py - extract reads based on coverage; see `khmer-recipes #1 <https://github.com/dib-lab/khmer-recipes/tree/master/001-extract-reads-by-coverage>`__

To keep, document, and build recipes for:

* `make-coverage.py - RPKM calculation script
* abundance-hist-by-position.py - look at abundance of k-mers by position within read; use with fasta-to-abundance-hist.py
* assemstats3.py - print out assembly statistics
* build-sparse-graph.py - code for building a sparse graph (by Camille Scott)
* calc-best-assembly.py - calculate the "best assembly" - used in metagenome protocol
* collect-variants.py - used in a `gist <https://gist.github.com/ctb/6eaef7971ea429ab348d>`__
* extract-single-partition.py - extract all the sequences that belong to a specific partition, from a file with multiple partitions
* fasta-to-abundance-hist.py - generate abundance of k-mers by position within reads; use with abundance-hist-by-position.py
* filter-below-abund.py - like filter-abund, but trim off high-abundance k-mers
* filter-median-and-pct.py - see blog post on Trinity in silico norm (http://ivory.idyll.org/blog/trinity-in-silico-normalize.html)
* filter-median.py - see blog post on Trinity in silico norm (http://ivory.idyll.org/blog/trinity-in-silico-normalize.html)
* find-high-abund-kmers.py - extract high-abundance k-mers into a list
* graph-size.py - filter reads based on size of connected graph
* hi-lo-abundance-by-position.py - look at high and low-abundance k-mers by position within read
* memusg - memory usage analysis
* multi-rename.py - rename sequences from multiple files with a common prefix
* normalize-by-median-pct.py - see blog post on Trinity in silico norm (http://ivory.idyll.org/blog/trinity-in-silico-normalize.html)
* print-stoptags.py - print out the stoptag k-mers
* print-tagset.py - print out the tagset k-mers
* renumber-partitions.py - systematically renumber partitions
* shuffle-reverse-rotary.py - FASTA file shuffler for larger FASTA files
* split-fasta.py - break a FASTA file up into smaller chunks
* stoptag-abundance-hist.py - print out abundance histogram of stoptags
* stoptags-by-position.py - print out where stoptags tend to occur
* strip-partition.py - clear off partition information
* subset-report.py - report stats on pmap files
* sweep-files.py - various ways to extract reads based on k-mer overlap
* sweep-out-reads-with-contigs.py - various ways to extract reads based on k-mer overlap
* sweep-reads.py - various ways to extract reads based on k-mer overlap
* sweep-reads2.py - various ways to extract reads based on k-mer overlap
* sweep-reads3.py - various ways to extract reads based on k-mer overlap
* write-trimmomatic.py - used to build Trimmomatic command lines in `khmer-protocols <http://khmer-protocols.readthedocs.org/en/latest/>`__

Good ideas to rewrite using newer tools/approaches:

* assembly-diff.py - find sequences that differ between two assemblies
* assembly-diff-2.py - find subsequences that differ between two assemblies
* bloom-count.py - count # of unique k-mers; should be reimplemented with HyperLogLog, Renamed from bloom_count.py in commit 4788c31
* bloom-count-intersection.py - look at unique and disjoint #s of k-mers, Renamed from bloom_count_intersection.py in commit 4788c31.
* split-sequences-by-length.py - break up short reads by length

----

Present in commit d295bc847 but removed thereafter:

* `combine-pe.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/combine-pe.py>`__ - combine partitions based on shared PE reads.
* `compare-partitions.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/compare-partitions.py>`__ - compare read membership in partitions.
* `count-within-radius.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/count-within-radius.py>`__ - calculating graph density by position with seq
* `degree-by-position.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/degree-by-position.py>`__ - calculating graph degree by position in seq
* `dn-identify-errors.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/dn-identify-errors.py>`__ - prototype script to identify errors in reads based on diginorm principles
* `ec.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/ec.py>`__ - new error correction foo
* `error-correct-pass2.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/error-correct-pass2.py>`__ - new error correction foo
* `find-unpart.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/find-unpart.py>`__ - something to do with finding unpartitioned sequences
* `normalize-by-align.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/normalize-by-align.py>`__  - new error correction foo
* `read_aligner.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/read_aligner.py>`__ - new error correction foo
* `shuffle-fasta.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/shuffle-fasta.py>`__ - FASTA file shuffler for small FASTA files
* `to-casava-1.8-fastq.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/to-casava-1.8-fastq.py>`__ - convert reads to different Casava format
* `uniqify-sequences.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/uniqify-sequences.py>`__ - print out paths that are unique in the graph
* `write-interleave.py <https://github.com/dib-lab/khmer/blob/d295bc8477022e8c34649f131a2abe333a891d3d/sandbox/write-interleave.py>`__ - is this used by any protocol etc?

Present in commit 691b0b3ae but removed thereafter:

* `annotate-with-median-count.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/annotate-with-median-count.py>`__ - replaced by count-median.py
* `assemble-individual-partitions.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/assemble-individual-partitions.py>`__ - better done with parallel
* `assemstats.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/assemstats.py>`__ - statistics gathering; see assemstats3.
* `assemstats2.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/assemstats2.py>`__ - statistics gathering; see assemstats3.
* `abund-ablate-reads.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/abund-ablate-reads.py>`__ - trim reads of high abundance k-mers.
* `bench-graphsize-orig.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/bench-graphsize-orig.py>`__ - benchmarking script for graphsize elimination
* `bench-graphsize-th.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/bench-graphsize-th.py>`__ - benchmarking script for graphsize elimination
* `bin-reads-by-abundance.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/bin-reads-by-abundance.py>`__ - see slice-reads-by-coverage.py
* `bowtie-parser.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/bowtie-parser.py>`__ - parse bowtie map file
* `calc-degree.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/calc-degree.py>`__ - various k-mer statistics
* `calc-kmer-partition-counts.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/calc-kmer-partition-counts.py>`__ - various k-mer statistics
* `calc-kmer-read-abunds.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/calc-kmer-read-abunds.py>`__ - various k-mer statistics
* `calc-kmer-read-stats.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/calc-kmer-read-stats.py>`__ - various k-mer statistics
* `calc-kmer-to-partition-ratio.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/calc-kmer-to-partition-ratio.py>`__ - various k-mer statistics
* `calc-sequence-entropy.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/calc-sequence-entropy.py>`__ - calculate per-sequence entropy
* `choose-largest-assembly.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/choose-largest-assembly.py>`__ - see calc-best-assembly.py
* `consume-and-traverse.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/consume-and-traverse.py>`__ - replaced by load-graph.py
* `contig-coverage.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/contig-coverage.py>`__ - calculate coverage of contigs by k-mers
* `count-circum-by-position.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/count-circum-by-position.py>`__ - k-mer graph statistics by position within read
* `count-density-by-position.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/count-density-by-position.py>`__ - k-mer graph stats by position within read
* `count-distance-to-volume.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/count-distance-to-volume.py>`__ - k-mer stats from graph
* `count-median-abund-by-partition.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/count-median-abund-by-partition.py>`__ - count median k-mer abundance by partition;
* `count-shared-kmers-btw-assemblies.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/count-shared-kmers-btw-assemblies.py>`__ - count shared k-mers between assemblies;
* `ctb-iterative-bench-2-old.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/ctb-iterative-bench-2-old.py>`__ - old benchmarking code
* `ctb-iterative-bench.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/ctb-iterative-bench.py>`__ - old benchmarking code
* `discard-high-abund.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/discard-high-abund.py>`__ - discard reads by coverage; see slice-reads-by-coverage.py
* `discard-pre-high-abund.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/discard-pre-high-abund.py>`__ - discard reads by coverage; see slice-reads-by-coverage.py
* `do-intertable-part.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-intertable-part.py>`__ - unused partitioning method
* `do-partition-2.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-partition-2.py>`__ - replaced by scripts/do-partition.py
* `do-partition-stop.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-partition-stop.py>`__ - replaced by scripts/do-partition.py
* `do-partition.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-partition.py>`__ - moved to scripts/
* `do-subset-merge.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-subset-merge.py>`__ - replaced by scripts/merge-partitions.py
* `do-th-subset-calc.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-th-subset-calc.py>`__ - unused benchmarking scripts
* `do-th-subset-load.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-th-subset-load.py>`__ - unused benchmarking scripts
* `do-th-subset-save.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/do-th-subset-save.py>`__ - unused benchmarking scripts
* `extract-surrender.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/extract-surrender.py>`__ - no longer used partitioning feature
* `extract-with-median-count.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/extract-with-median-count.py>`__ - see slice-reads-by-coverage.py
* `fasta-to-fastq.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/fasta-to-fastq.py>`__ - just a bad idea
* `filter-above-median.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-above-median.py>`__ - replaced by filter-below-abund.py
* `filter-abund-output-by-length.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-abund-output-by-length.py>`__ - replaced by filter-abund/filter-below-abund
* `filter-area.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-area.py>`__ - trim highly connected k-mers
* `filter-degree.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-degree.py>`__ - trim highly connected k-mers
* `filter-density-explosion.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-density-explosion.py>`__ - trim highly connected k-mers
* `filter-if-present.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-if-present.py>`__ - replaced by filter-abund and others
* `filter-max255.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-max255.py>`__ - remove reads w/high-abundance k-mers.
* `filter-min2-multi.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-min2-multi.py>`__ - remove reads w/low-abundance k-mers
* `filter-sodd.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-sodd.py>`__ - no longer used partitioning feature
* `filter-subsets-by-partsize.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/filter-subsets-by-partsize.py>`__ - deprecated way to filter out partitions by size
* `get-occupancy.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/get-occupancy.py>`__ - utility script no longer needed
* `get-occupancy2.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/get-occupancy2.py>`__ - utility script no longer needed
* `graph-partition-separate.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/graph-partition-separate.py>`__ - deprecated graph partitioning stuff
* `graph-size-circum-trim.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/graph-size-circum-trim.py>`__ - experimental mods to graph-size.py
* `graph-size-degree-trim.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/graph-size-degree-trim.py>`__ - experimental mods to graph-size.py
* `graph-size-py.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/graph-size-py.py>`__ - experimental mods to graph-size.py
* `join_pe.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/join_pe.py>`__ - silly attempts to deal with PE interleaving?
* `keep-stoptags.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/keep-stoptags.py>`__ - trim at stoptags
* `label-pairs.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/label-pairs.py>`__ - deprecated PE fixing script
* `length-dist.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/length-dist.py>`__ - deprecated length distribution calc script
* `load-ht-and-tags.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/load-ht-and-tags.py>`__ - load and examine hashtable & tags
* `multi-abyss.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/multi-abyss.py>`__ - better done with parallel
* `make-coverage-by-position-for-node.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/make-coverage-by-position-for-node.py>`__ - deprecated coverage calculation
* `make-coverage-histogram.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/make-coverage-histogram.py>`__ - build coverage histograms
* `make-random.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/make-random.py>`__ - make random DNA; see dbg-graph-null project.
* `make-read-stats.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/make-read-stats.py>`__ - see readstats.py
* `multi-stats.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/multi-stats.py>`__ - see readstats.py
* `multi-velvet.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/multi-velvet.py>`__ - better done with parallel
* `normalize-by-min.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/normalize-by-min.py>`__ - normalize by min k-mer abundance in seq; just a bad idea
* `occupy.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/occupy.py>`__ - no longer needed utility script
* `parse-bowtie-pe.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/parse-bowtie-pe.py>`__ - no longer needed utility script
* `parse-stats.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/parse-stats.py>`__ - partition stats
* `partition-by-contig.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/partition-by-contig.py>`__ - various approaches to partitioning
* `partition-by-contig2.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/partition-by-contig2.py>`__ - various approaches to partitioning
* `partition-size-dist-running.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/partition-size-dist-running.py>`__ - various approaches to partitioning
* `partition-size-dist.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/partition-size-dist.py>`__ - various approaches to partitioning
* `path-compare-to-vectors.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/path-compare-to-vectors.py>`__ - ??
* `print-exact-abund-kmer.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/print-exact-abund-kmer.py>`__ - ??
* `print-high-density-kmers.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/print-high-density-kmers.py>`__ - display high abundance k-mers
* `quality-trim-pe.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/quality-trim-pe.py>`__ - no longer needed utility script
* `quality-trim.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/quality-trim.py>`__ - no longer needed utility script
* `reformat.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/reformat.py>`__ - FASTA sequence description line reformatter for partitioned files
* `remove-N.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/remove-N.py>`__ - eliminate sequences that have Ns in them
* `softmask-high-abund.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/softmask-high-abund.py>`__ - softmask high abundance sequences (convert ACGT to acgt)
* `split-fasta-on-circum.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-circum.py>`__ - various ways of breaking sequences on graph properties
* `split-fasta-on-circum2.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-circum2.py>`__ - various ways of breaking sequences on graph properties
* `split-fasta-on-circum3.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-circum3.py>`__ - various ways of breaking sequences on graph properties
* `split-fasta-on-circum4.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-circum4.py>`__ - various ways of breaking sequences on graph properties
* `split-fasta-on-degree-th.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-degree-th.py>`__ - various ways of breaking sequences on graph properties
* `split-fasta-on-degree.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-degree.py>`__ - various ways of breaking sequences on graph properties
* `split-fasta-on-density.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-fasta-on-density.py>`__ - various ways of breaking sequences on graph properties
* `split-N.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-N.py>`__ - truncate sequences on N
* `split-reads-on-median-diff.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/split-reads-on-median-diff.py>`__ - various ways of breaking sequences on graph properties
* `summarize.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/summarize.py>`__ - sequence stats calculator
* `sweep_perf.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/sweep_perf.py>`__ - benchmarking tool
* `test_scripts.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/test_scripts.py>`__ - old test file
* `traverse-contigs.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/traverse-contigs.py>`__ - deprecated graph traversal stuff
* `traverse-from-reads.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/traverse-from-reads.py>`__ - deprecated graph traversal stuff
* `validate-partitioning.py <https://github.com/dib-lab/khmer/tree/691b0b3aefe83e9e8f5f2b80f5f9516664a4654a/sandbox/validate-partitioning.py>`__ - unneeded test
