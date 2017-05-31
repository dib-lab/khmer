# khmer v2.0 release notes

This is the v2.0 release of khmer and the first from our new lab at the
University of California, Davis. It features Python 3 compatibility, streaming
I/O from Unix Pipes, mixed-pair sequence file format support, and a new
parameter to simplify memory usage. We also have a software paper in-press
describing the project and the citation reminders have been updated to reflect
that.

Overall there are an additional 2,380 lines of Python code (mostly tests) and
283 less lines of C++ (despite adding features). This release is the product of
over 1,000 commits to the codebase since v1.4.

Documentation is at https://khmer.readthedocs.org/en/v2.0/

## New items of note:

New behavior
------------

### Streaming I/O from Unix Pipes

All scripts now accept input from named (like `/dev/stdin`, or that
created using `<( list )` process substituion) and unnamed pipes (like
output piped in from another program with `|`). The STDIN stream can
also be specified using a single dash: `-`. #1186 @mr-c #1042 #763 @SherineAwad
\#1085 @ctb

### New parameter for memory usage, and/or tablesize/number of table parameters.

There is now a `-M`/`--max-memory-usage` parameter that sets the number of
tables (`-N`/`--n\_tables`) and tablesize (`-x`/`--max-tablesize`) parameters
automatically to match the desired memory usage. #1106 #621 #1126 #390 #1117
\#1055 #1050 #1214 #1179 #1133 #1145 @ctb @qingpeng @bocajnotnef  

### Digital normalization script now supports mixed paired and unpaired read input

`normalize-by-median.py` now supports mixed paired and unpaired (or
"broken-paired") input. Behavior can be forced to either treat all reads
as singletons or to require all reads be properly paired
using `--force\_single` or `--paired`, respectively. If `--paired` is
set, `--unpaired-reads` can be used to include a file of unpaired reads. The
unpaired reads will be examined after all of the other sequence files.
`normalize-by-median.py` now has a `--quiet` option to reduce the amount of
output. #1200 @bocajnotnef

### Mixed-pair sequence file format support

`split-paired-reads.py` `--output-orphaned`/`-0` has been added to allow for
orphaned reads and give them a file to be sorted into. #847 #1164 @ctb

### Scripts now output columnar data in CSV format by default

All scripts that output any kind of columnar data now do so in CSV
format, with headers. Previously this had to be enabled with `--csv`.
(Affects `abundance-dist-single.py`, `abundance-dist.py`, `count-median.py`,
and `count-overlap.py`.) `normalize-by-median.py --report` also now outputs
in CSV format. #1011 #1180 @ctb

### Reservoir sampling script extracts paired reads by default

`sample-reads-randomly.py` now retains pairs in the output, by default.
This can be overridden to match previous behavior with `--force\_single`.

### Most input and output files can be compressed

We support gzip and bzip2 input and output file compression everywhere that it
makes sense #505 #747 @bocajnotnef

New scripts
-----------

### Estimate number of unique kmers

`unique-kmers.py` estimates the k-mer cardinality of a dataset using the
HyperLogLog probabilistic data structure. This allows very low memory
consumption, which can be configured through an expected error rate.
Even with low error rate (and higher memory consumption), it is still
much more efficient than exact counting and alternative methods. It
supports multicore processing (using OpenMP) and streaming, and so can
be used in conjunction with other scripts (like `normalize-by-median.py`
and `filter-abund.py`). This script is the work of @luizirber and the subject
of a paper in draft. #390 #1239 #1252 #1053 #1072 #1145 #1176 #1207 #1204 #1245

Incompatible changes
--------------------

### New datastructure and script names

For clarity the Count-Min Sketch based data structure previously known
as "counting\_hash" or "counting\_table" and variations of these is now
known as `countgraph`. Likewise with the Bloom Filter based data
structure previously known at "hashbits", "presence\_table" and
variations of these is now known as `nodegraph`. Many options relating
to `table` have been changed to `graph`. #1112 #1209 @mr-c

### Binary file formats have changed

All binary khmer formats (presence tables, counting tables, tag sets,
stop tags, and partition subsets) have changed. Files are now pre-pended
with the string `OXLI` to indicate that they are from this project. #519 #1031
@mr-c #1159 @luizirber

Files of the above types made in previous versions of khmer are not
compatible with v2.0; the reverse is also true.

In addition to the `OXLI` string, the Nodegraph and Countgraph file
format now includes the number of occupied bins. See
http://khmer.readthedocs.org/en/v2.0/dev/binary-file-formats for details. #1093
@ctb @mr-c #1101 #1103 @kdmurray91

### load-graph.py no longer appends .pt to the specified filename

Previously, `load-graph.py` appended a `.pt` extension to the specified
output filename and partition-graph.py appended a `.pt` to the given
input filename. Now, `load-graph.py` writes to the specified output
filename and `partition-graph.py` does not append a `.pt` to the given
input filename. #1189 #747 @bocajnotnef

### Some reporting options have been turned always on

The total number of unique k-mers will always be reported every time a
new countgraph is made. The `--report-total-kmers` option has been
removed from `abundance-dist-single.py`, `filter-abund-single.py`, and
`normalize-by-median.py` to reflect this. Likewise with `--write-fp-rate`
for `load-into-counting.py` and `load-graph.py`; the false positive rate
will always be written to the `.info` files. #1097 #1180 @ctb

### An uncommon error recovery routine was removed

To simplify the codebase `--save-on-failure` and its helper option
`--dump-frequency` have been removed from `normalize-by-median.py`.

### Single file output option names have been normalized

`--out` is now `--output` for both `normalize-by-median.py` and 
`trim-low-abund.py`. #1188 #1164 @ctb

### Miscellaneous changes

The common option `--min-tablesize` was renamed to `--max-tablesize` to reflect
this more desirable behavior.

In conjuction with the new `split-paired-reads.py` `--output-orphaned`
option, the option `--force-paired`/`-p` has been eliminated.

As CSV format is now the default, the `--csv` option has been removed.

### Removed script

[count-overlap.py](http://khmer.readthedocs.org/en/v1.4.1/user/scripts.html#count-overlap-py)
has been removed.

## Notable bugs fixed/issues closed:

When `normalize-by-median.py` decides to keep both parts of a pair of reads it
was only adding the k-mers & counts from one to the countgraph. #1000 #1010
@drtamermansour @bocajnotnef

The partition map file format was not robust to truncation and would hang
waiting for more data. #437 #1037 #1048 @ctb

`extract-paired-reads.py` and `split-paired-reads.py` no longer create default
files when the user supplies filename(s). #1005 #1132 @kdmurray91

## Additional fixes/features

`find-knots.py` was missing a `--force` option and unit tests. #358 #1078 @ctb
The check for excessively high false-positive rate has also received a
`--force` option #1168 @bocajnotnef

A bug leading to an infinite loop with large gzipped countgraphs was found
\#1038 #1043 @kdmurray91

All scripts that create nodegraphs or countgraphs report the total number of
unique k-mers. #491 #609 #429 @mr-c

Read pairs from SRA are fully supported. Reported by @macmanes in #1027, fixed
by @kdmurray91 @SherineAwad in #1173 #1088

### Of interest to users:

Added `Hashtable::get_kmers()`, `get_kmer_hashes()`, and `get_kmer_counts()`
with corresponding CPython functions. #1047 #1049 @ctb 

The `DEFAULT_DESIRED_COVERAGE` for `normalize-by-median.py` is now 20. #1073
\#1081 @ctb

FIFOs are no longer seen as empty. #1147 #1163 @bocajnotnef

When the k-size is requested to be larger than 32 (which is unsupported) a
helpful error message is reported. #1094 #1050 @ctb

We try to report more helpfully during errors, such as suggesting the `--force`
option when outputs files already exist. #1162 #1170 @bocajnotnef

There is a paper related to `trim-low-abund.py`: "Crossing the streams: a
framework for streaming analysis of short DNA sequencing reads" and it has been
added to the CITATION file and program output. #1180 #1130 @ctb

We have dropped support for Python 2.6 #1009 #1180 @ctb

Our user documentation got a bit out of date and has been updated. #1156 #1247
@bocajnotnef @mr-c #1104 @kdmurray91 #1267 @ctb Links to lists of publications
that use khmer have been added #1063 #1222 @mr-c The help text from the scripts
has also had a thorough cleanup for formatting. #1268 @mr-c

`fastq-to-fasta.py`'s `--n_keep` option has incorrect help text. We now point
out that all reads with Ns will be dropped by default unless this option is
supplied. #657 #814 #1208 @ACharbonneau @bocajnotnef

We've updated the URL to the '88m-reads.fa.gz' file. #1242 #1269 @mr-c

@camillescott designed and implemented an optimization for
`normalize-by-median.py` #862

`abundance-dist.py` can now be used without counts over 255 with
`--no-bigcount`. #1067 #909 @drtamermansour @bocajnotnef Its input file
requirement can no longer be overridden #1201 #1202 @bocajnotnef

khmer v2.0 will be released as a package for the Debian GNU/Linux operating
system. Big thanks to @kdmurray91 for his assistance. #1148 #1240 The C++
library, now named liboxli, will have its own package as well.

`sandbox/multi-rename.py` now wraps long FASTA sequences at 80 columns. #450
\#1136 @SherineAwad

### Of interest to developers:

The khmer project is now a Python 3 codebase with backwards compatibility to
Python 2.7. Huge credit to @luizirber #978 #922 #1045 #1066 #1089 #1157 #1191
\#1108 Many developer impacting changes including the file
`khmer/\_khmermodule.cc` is now `khmer/\_khmer.cc`. #169 #904

@camillescott did an extensive refactor of the C++ graph traversal code which
removed a considerable amount of redundant code and will be very useful for
future work. #1231 #1080

We now use some and allow all C++11 features in the codebase. #598 #1122 @mr-c 

`normalize-by-median.py` was extensively refactored. #1006 #1010 #1057 #1039
\#1135 #1182 @bocajnotnef @ctb @camillescott

The CPython glue was refactored so that CountingHash and Hashbits inherit from
Hashtable. #1044 @ctb

The tests no longer stop on the first failed test. #1124 #1134 @ctb and some
noisy tests were silenced #1125 #1137 @bocajnotnef

The `check_space()` calls were cleaned up. #1167 #1166 #1170 #993 

Developer docs have been expanded #737 #1184 @bocajnotnef #1083 #1282 @ctb
@mr-c #1269

A lot of code was deleted: TRACE related code in #274 #1180 @ctb 
`hashtable_collect_high_abundance_kmers` in #1142 #1044 @ctb `lib/ht-diff.cc`,
`lib/test-HashTables.cc`, `lib/test-Parser.cc` #1144, @mr-c `bink.ipynb`,
`lib/graphtest.cc`, `lib/primes.hh` #1289 @mr-c

@bocajnotnef deleted more unused code and added new tests elsewhere to increase
testing coverage in #1236. @mr-c had his own go in #1279

cppcheck installation for OSX has been documented #777 #952 #945 @elmbeech

ccache and git-merge-changelog has been documented for Linux users #610 #1122
\#614 @mr-c

The graphalign parameters can be saved/loaded from disk. In addition the
`align_forward` method has been introduced. #755 #750 @mr-c @ctb

`labelhash` is now known as `graphlabels` #1032 #1209 @mr-c It is also now a
'friend' of Hashtable and one can make either a nodegraph or countgraph
version. These graphlabels can now be saved & loaded from disk. #1021 @ctb

Spelling is hard; we've added instructions on how to run codespell to the
developer docs. #890 #1203 @bocajnotnef

A redundant and contradictory named test has been removed. Reported by @jgluck
in #662 fixed by @bocajnotnef in #1220 @SherineAwad contributed some additional
tests #809 #615.

The new oxli command, while disabled in the v2.0 release, has been added to all
the QA makefile targets as we continue to refactor the codebase. #1199 #1218
@bocajnotnef

The CPython code was audited to ensure that all possible C++ exceptions were
caught and dealt with. The exception hierarchy was also simplified #1016 #1015
\#1017 #1151 @kdmurray91 @mr-c

`get_kadian_count` has been removed. #1034 #1194 @ctb

We use argparse's `metavar`s to aid with autogenerated documentation for the
scripts. This has been documented in the dev docs. #620 #1222 @mr-c

Sometimes one makes a lot of commits while refining a feature or pull request.
We've documented a field-tested way to turn a pile of commits into a single
commit without the pain of `git rebase`. #1013 #660 #1222 @mr-c

We use Coverity to test for various issues with our C++ code. The Makefile
target has been updated for changes on their side. #1007 #1222 @mr-c

There is a new `update()` function to merge two nodegraphs of the same size and
ksize. #1051 @ctb

Despite the checklist, formatting errors still occur. We must be vigilant!
\#1075 @luizirber

There is a new `filter_on_median` function. #862 #1077 @camillescott

There are new scripts in the `sandbox/` which output k-mer counts:
sandbox/{count-kmers.py,count-kmers-single.py}. #983 @ctb

A large effort to make the codebase 'pylint clean' has begun with #1175
@bocajnotnef Likewise the cpychecker tool was re-run on the CPython code and
issues found there were addressed #1196 @mr-c

As repeatedly promised, we've updated our list of contributors to include
everyone with a commit in git. #1023 @mr-c

`thread_utils.is_pair()` has been dropped in favor of `utils.check_is_pair()`
\#1284 @mr-c

The Doxygen produced documentation is improving. The location of included
headers is now autodetected for Doxygen and cppcheck.

## Known issues:

``load-graph.py`` in multithreaded mode will find slightly different number of
unique kmers. This is being investigated in #1248

## Contributors

@ctb, @bocajnotnef, @mr-c, @luizirber, @kdmurray91, @SherineAwad,
@camillescott, ‡@ACharbonneau

‡ Indicates new contributors

## Issue reporters

@jgluck, @ACharbonneau, @macmanes
