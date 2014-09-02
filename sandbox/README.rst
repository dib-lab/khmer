Sandbox scripts
===============

To keep:

* abundance-hist-by-position.py - look at abundance of k-mers by position within read; use with fasta-to-abundance-hist.py
* assemstats3.py - print out assembly statistics
* calc-best-assembly.py - calculate the "best assembly" - used in metagenome protocol
* calc-median-distribution.py - plot coverage distribution; see `khmer-recipes #1 <https://github.com/ged-lab/khmer-recipes/tree/master/001-extract-reads-by-coverage>`__
* combine-pe.py - combine partitions based on shared PE reads.
* compare-partitions.py
* dn-identify-errors.py - prototype script to identify errors in reads based on diginorm principles
* extract-single-partition.py - extract all the sequences that belong to a specific partition, from a file with multiple partitions
* fasta-to-abundance-hist.py - generate abundance of k-mers by position within reads; use with abundance-hist-by-position.py
* filter-median-and-pct.py - see blog post on Trinity in silico norm (http://ivory.idyll.org/blog/trinity-in-silico-normalize.html)
* filter-median.py - see blog post on Trinity in silico norm (http://ivory.idyll.org/blog/trinity-in-silico-normalize.html)
* graph-size.py - filter reads based on size of connected graph
* multi-rename.py - rename sequences from multiple files with a common prefix
* normalize-by-median-pct.py - see blog post on Trinity in silico norm (http://ivory.idyll.org/blog/trinity-in-silico-normalize.html)
* readstats.py - print out read statistics
* saturate-by-median.py - calculate collector's curve on shotgun sequencing; see `khmer-recipes #4 <https://github.com/ged-lab/khmer-recipes/tree/master/004-estimate-sequencing-saturation>`__
* trim-low-abund.py - streaming version of filter-abund.

Good ideas to rewrite using newer tools/approaches::

* assembly-diff.py - find sequences that differ between two assemblies
* assembly-diff-2.py - find subsequences that differ between two assemblies
* bloom_count.py - count # of unique k-mers; should be reimplemented with HyperLogLog
* bloom_count_intersection.py - look at unique and disjoint #s of k-mers

To examine:

* build-sparse-graph.py - code for building a sparse graph
* count-within-radius.py
* degree-by-position.py
* do-subset-merge.py
* ec.py
* error-correct-pass2.py
* normalize-by-align.py
* normalize-by-min.py

Of unknown utility:

* filter-abund-output-by-length.py
* filter-area.py
* filter-below-abund.py
* filter-degree.py
* filter-density-explosion.py
* filter-if-present.py
* filter-subsets-by-partsize.py
* find-high-abund-kmers.py
* find-unpart.py
* graph-partition-separate.py
* hi-lo-abundance-by-position.py
* keep-stoptags.py
* label-pairs.py
* length-dist.py
* make-coverage-by-position-for-node.py
* make-coverage-histogram.py
* make-coverage.py
* partition-by-contig.py
* partition-by-contig2.py
* partition-size-dist-running.py
* partition-size-dist.py
* path-compare-to-vectors.py
* print-exact-abund-kmer.py
* print-high-density-kmers.py
* print-stoptags.py
* print-tagset.py
* read_aligner.py
* reformat.py
* remove-N.py
* renumber-partitions.py
* shuffle-fasta.py
* shuffle-reverse-rotary.py
* softmask-high-abund.py
* split-N.py
* split-fasta-on-circum.py
* split-fasta-on-circum2.py
* split-fasta-on-circum3.py
* split-fasta-on-circum4.py
* split-fasta-on-degree-th.py
* split-fasta-on-degree.py
* split-fasta-on-density.py
* split-fasta.py
* split-reads-on-median-diff.py
* split-sequences-by-length.py
* stoptag-abundance-hist.py
* stoptags-by-position.py
* strip-partition.py
* subset-report.py
* summarize.py
* sweep-files.py
* sweep-out-reads-with-contigs.py
* sweep-reads.py
* sweep-reads2.py
* sweep-reads3.py
* sweep_perf.py
* test_scripts.py
* to-casava-1.8-fastq.py
* traverse-contigs.py
* traverse-from-reads.py
* uniqify-sequences.py
* validate-partitioning.py
* write-interleave.py
* write-trimmomatic.py

Present in commit 691b0b3ae but removed thereafter:

* annotate-with-median-count.py - replaced by count-median.py
* assemble-individual-partitions.py - better done with parallel
* assemstats.py - statistics gathering; see assemstats3.
* assemstats2.py - statistics gathering; see assemstats3.
* abund-ablate-reads.py - trim reads of high abundance k-mers.
* bench-graphsize-orig.py - benchmarking script for graphsize elimination
* bench-graphsize-th.py - benchmarking script for graphsize elimination
* bin-reads-by-abundance.py - see slice-reads-by-coverage.py
* bowtie-parser.py - parse bowtie map file
* calc-degree.py - various k-mer statistics
* calc-kmer-partition-counts.py - various k-mer statistics
* calc-kmer-read-abunds.py - various k-mer statistics
* calc-kmer-read-stats.py - various k-mer statistics
* calc-kmer-to-partition-ratio.py - various k-mer statistics
* calc-sequence-entropy.py - calculate per-sequence entropy
* choose-largest-assembly.py - see calc-best-assembly.py
* consume-and-traverse.py - replaced by load-graph.py
* contig-coverage.py - calculate coverage of contigs by k-mers
* count-circum-by-position.py - k-mer graph statistics by position within read
* count-density-by-position.py - k-mer graph stats by position within read
* count-distance-to-volume.py - k-mer stats from graph
* count-median-abund-by-partition.py - count median k-mer abundance by partition;
* count-shared-kmers-btw-assemblies.py - count shared k-mers between assemblies;
* ctb-iterative-bench-2-old.py - old benchmarking code
* ctb-iterative-bench.py - old benchmarking code
* discard-high-abund.py - discard reads by coverage; see slice-reads-by-coverage.py
* discard-pre-high-abund.py - discard reads by coverage; see slice-reads-by-coverage.py
* do-intertable-part.py - unused partitioning method
* do-partition-2.py - replaced by scripts/do-partition.py
* do-partition-stop.py - replaced by scripts/do-partition.py
* do-partition.py - moved to scripts/
* do-th-subset-calc.py - unused benchmarking scripts
* do-th-subset-load.py - unused benchmarking scripts
* do-th-subset-save.py - unused benchmarking scripts
* extract-surrender.py - no longer used partitioning feature
* extract-with-median-count.py - see slice-reads-by-coverage.py
* fasta-to-fastq.py - just a bad idea
* filter-above-median.py - replaced by filter-below-abund.py
* filter-max255.py - remove reads w/high-abundance k-mers.
* filter-min2-multi.py - remove reads w/low-abundance k-mers
* filter-sodd.py - no longer used partitioning feature
* get-occupancy.py - utility script no longer needed
* get-occupancy2.py - utility script no longer needed
* graph-size-circum-trim.py - experimental mods to graph-size.py
* graph-size-degree-trim.py - experimental mods to graph-size.py
* graph-size-py.py - experimental mods to graph-size.py
* join_pe.py - silly attempts to deal with PE interleaving?
* load-ht-and-tags.py - load and examine hashtable & tags
* multi-abyss.py - better done with parallel
* make-random.py - make random DNA; see dbg-graph-null project.
* make-read-stats.py - see readstats.py
* multi-stats.py - see readstats.py
* multi-velvet.py - better done with parallel
* occupy.py - no longer needed utility script
* parse-bowtie-pe.py - no longer needed utility script
* parse-stats.py - partition stats
* quality-trim-pe.py - no longer needed utility script
* quality-trim.py - no longer needed utility script
