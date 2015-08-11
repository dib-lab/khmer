# khmer v1.1 release notes

This is v1.1, a minor version release; this version adds several new scripts.

Docs at: https://khmer.readthedocs.org/en/v1.1/

Release notes w/links: https://github.com/dib-lab/khmer/releases/tag/v1.1

## New items of note:

* removed unnecessary files from PyPI package; distribution is now under 2 MB (#419) @mr-c
* tests are now distributed with package and can be run after 'pip install' (#451) @mr-c
* complain properly on file read failures (#333) @ctb
* Sequence loading scripts will now report total numbers of k-mers if given --report_total_kmers (#491/#429) @mr-c
* added metagenome protocol to acceptance testing (#472) @SherineAwad @ctb

## Notable bugs fixed/issues closed:

* removed sandbox/load-into-hashbits.py (superseded by scripts/load-graph.py --no-tagset) (#480, @wrightmhw)
* promoted extract-long-sequences.py to scripts (#461, @wrightmhw)
* promoted fastq-to-fasta.py to scripts (#436, @wrightmhw)
* remove incorrect filesystem space check from abundance-dist.py (#452, @chuckpr)
* when counting hash writes fail, produce error message (#411, @znruss)
* removed a number of memory leaks found by Coverity and valgrind (#451, @mr-c)
* updated reservoir sampling to produce multiple subsamples with -S (#197, @ctb)
* fixed pip2, python2 issues (#428 and #485, @accaldwell @mr-c)
* removed untested/unused code and scripts (#438, @mr-c)

## Known issues:

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/dib-lab/khmer/issues/266

Some FASTQ files confuse our parser when running with more than one thread.
For example, while using load-into-counting.py. If you experience this then
add "--threads=1" to your command line. This issue is being tracked in
https://github.com/dib-lab/khmer/issues/249

If your k-mer table is truncated on write, an error may not be reported; this
is being tracked in https://github.com/dib-lab/khmer/issues/443.
However, khmer will now (correctly) fail when trying to read a truncated file
(See #333).

Paired-end reads from Casava 1.8 currently require renaming for use in
normalize-by-median and abund-filter when used in paired mode. The
integration of a fix for this is being tracked in https://github.com/dib-lab/khmer/issues/23

Some scripts only output FASTA even if given a FASTQ file. This issue
is being tracked in https://github.com/dib-lab/khmer/issues/46

A user reported that abundance-dist-single.py fails with small files and many
threads. This issue is being tracked in https://github.com/dib-lab/khmer/issues/75

## Contributors

@mr-c, @ctb, @camillescott, @wrightmhw, @chuckpr, @luizirber, @accaldwell,
@znruss
