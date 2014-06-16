This is a bug fix release.

Docs at: https://khmer.readthedocs.org/en/v1.1/

Release notes w/links: https://github.com/ged-lab/khmer/releases/tag/v1.1

## New items of note:

removed unnecessary files from PyPI package; distribution is now under 2 MB (#419)
tests are now distributed with package and can be run after 'pip install' (#451)
complain properly on file read failures (#333)

## Notable bugs fixed/issues closed:

removed sandbox/load-into-hashbits.py (superseded by scripts/load-graph.py --no-tagset) (#480)
Sequence loading scripts will now report total numbers of k-mers if given --report_total_kmers (#491/#429)
promoted extract-long-sequences.py to scripts (#461)
promoted fastq-to-fasta.py to scripts (#436)
removed a number of memory leaks found by Coverity and valgrind (#451)
updated reservoir sampling to produce multiple subsamples with -S (#197)
fixed pip2, python2 issues (#428 and #485)
removed untested/unused code and scripts (#438)

## Known issues:

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/ged-lab/khmer/issues/266

Some FASTQ files confuse our parser when running with more than one thread.
For example, while using load-into-counting.py. If you experience this then
add "--threads=1" to your command line. This issue is being tracked in
https://github.com/ged-lab/khmer/issues/249

If your k-mer table (hashfile) gets truncated, perhaps from a full filesystem, then our
tools currently will get stuck. This is being tracked in https://github.com/ged-lab/khmer/issues/247 and https://github.com/ged-lab/khmer/issues/246

Paired-end reads from Casava 1.8 currently require renaming for use in
normalize-by-median and abund-filter when used in paired mode. The
integration of a fix for this is being tracked in https://github.com/ged-lab/khmer/issues/23

annotate-partitions.py only outputs FASTA even if given a FASTQ file. This
issue is being tracked in https://github.com/ged-lab/khmer/issues/46

A user reported that abundance-dist-single.py fails with small files and many
threads. This issue is being tracked in https://github.com/ged-lab/khmer/issues/75

## Contributors

@mr-c, @ctb, @camillescott, @wrightmhw, @chuckpr, @luizirber, @accaldwell,
@znruss
