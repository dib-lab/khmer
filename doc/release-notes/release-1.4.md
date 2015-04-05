# khmer v1.4 release notes

This is the v1.4 release of khmer featuring the results of our March coding sprint; the use of the new v0.8 release of screed (the library we use for pure Python reading of nucleotide sequence files); and the addition of @luizirber's HyperLogLog counter for quick cardinality estimation.

Documentation is at https://khmer.readthedocs.org/en/v1.4/

## New items of note:

Casava 1.8 read naming is now fully supported and in general the scripts no longer mangle read names. Side benefits: `split-paired-reads.py` will no longer drop reads with 'bad' names; `count-median.py` can generate output in CSV format. #759 #818 @ctb

Most scripts now support a "broken" interleaved paired-read format for FASTA/FASTQ nucleotide sequence files. [`trim-low-abund.py`](http://khmer.readthedocs.org/en/v1.4/user/scripts.html#trim-low-abund-py) has been promoted from the sandbox as well. #759 @ctb 

The script to transform an interleaved paired-read nucleotide sequence file into two now allows one to name the output files which can be useful in combination with named pipes for streaming processing #762 @ctb 

Streaming everywhere: thanks to screed v0.8 we now support streaming of almost all inputs and outputs. #830 @aditi9783 #812 @mr-c

Need a quick way to count total number of unique k-mers in very low memory? the `unique-kmers.py` in the sandbox uses a HyperLogLog counter to quickly (and with little memory) provide an estimate with a controllable error rate. #257 #738 #895 #902 @luizirber 

## Notable bugs fixed/issues closed:

Paired-end reads from Casava 1.8 no longer require renaming for use in normalize-by-median and abund-filter when used in paired mode #818 @ctb

Python version support clarified. We do not (yet) support Python 3.x #741 @mr-c

If a single output file mode is chosen for normalize-by-median.py we now default to overwriting the output. Appending the output is available by using the append redirection operator from the shell. #843 @drtamermansour 

Scripts that consume sequence data using C++ will now properly throw an error on truncated files. #897 @kdmurray91 

## Additional fixes/features

### Of interest to users:

Misc documentation updates #753 @PamelaM, #782 @bocajnotnef, #845 @alameldin, #804 @ctb, #870 @SchwarzEM

Installation instructions for Arch Linux have been added #723 @reedacartwright

The example script for the STAMPS database has been fixed to run correctly #781 @drtamermansour 

split-paired-reads.py: added -o option to allow specification of an output directory #752 @bede 

Fixed a string formatting glitch in `sample-reads-randomly.py` #773 @qingpeng 

CSV output also added to `abundance-dist.py`, `abundance-dist-single.py`, and `count-overlap.py` #831 #854 #855 @drtamermansour 

interleave-reads.py now prints the output filename nicely #827 @kdmurray91 

Cleaned up error for input file not existing #772 @jessicamizzi #851 @ctb 

Fixed error in `find-knots.py` #860 @TheOneHyer 

The help text for `load-into-counting.py` for the `--no-bigcounts`/`-b` flag has been clarified #857 @kdmurray91

### Of interest to developers:

Switched away from using `--user` install for developers #740 @mr-c @drtamermansour & #883 @standage 

Developers can now see a summary of important Makefile targets via `make help` #783 @standage 

The unused `khmer.load_pe` module has been removed #828 @kdmurray91 

Versioneer bug due to new screed release was squashed #835 @mr-c

A Python 2.6 and 2.7.2 specific bug was worked around #869 @kdmurray91 @ctb 

added functions hash_find_all_tags_list and hash_get_tags_and_positions to CountingHash objects #749 #765 @ctb

The `make diff-cover` and ChangeLog formatting requirements have been added to checklist #766 @mr-c 

A useful message is now presented if large tables fail to allocate enough memory #704 @mr-c

A checklist for developers adding new CPython types was added #727 @mr-c

Specific policies for sandbox/ and scripts/ content, and a process for adding new command line scripts into scripts/ have been added to the developer documentation #799 @ctb

Sandbox scripts update: corrected #! Python invocation #815 @Echelon9, executable bits, copyright headers,  no underscores in filenames #823 #826 #850 @alameldin several scripts deleted, docs + requirements updated #852 @ctb

Avoid running big-memory tests on OS X #819 @ctb

Unused callback code was removed #698 @mr-c

The CPython code was updated to use the new checklist and follow additional best practices #785 #842 @luizirber 

Added a read-only view of the raw counting tables #671 @camillescott #869 @kdmurray91 

Added a Python method for quickly getting the number of underlying tables in a counting or presence table #879 #880 @kdmurray91 

The C++ library can now be built separately for the brave and curious developer #788 @kdmurray91 

The ReadParser object now keeps track of the number of reads processed #877 @kdmurray91 

Documentation is now reproducible #886 @mr-c

Python future proofing: specify floor division #863 @mr-c

Miscellaneous spelling fixes; thanks codespell! #867 @mr-c

## Known issues:

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more memory than it was configured for. This is being investigated in https://github.com/ged-lab/khmer/issues/266

If your k-mer table is truncated on write, an error may not be reported; this is being tracked in https://github.com/ged-lab/khmer/issues/443. However, khmer will now (correctly) fail when trying to read a truncated file (See #333).

Some scripts only output FASTA even if given a FASTQ file. This issue is being tracked in https://github.com/ged-lab/khmer/issues/46

A user reported that abundance-dist-single.py fails with small files and many threads. This issue is being tracked in https://github.com/ged-lab/khmer/issues/75

## Contributors

@ctb, @kdmurray91, @drtamermansour, @luizirber, @mr-c, @jessicamizzi,
@standage, @bocanotnef, \*@alameldin, \*@aditi9783, \*@TheOneHyer, \*@bede,
\*@SchwarzEM, \*@reedacartwright, @Echelon9, @qingpeng, @SherineAwad, @PamelaM 

\* Indicates new contributors

## Issue reporters

@moorepants, @teshomem, @macmanes, @lexnederbragt, @r-gaia-cs

