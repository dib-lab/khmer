khmer v1.0 release notes
========================

582 changed files with 40,527 additions and 31,772 deletions.

The team has been hard at work since v0.8 to refine the codebase into a
stable product.

https://khmer.readthedocs.org/en/latest/

With the 1.0 release we are making a commitment to using Semantic
Versioning[0]: the version number will reflect the impact of the changes
between releases. New major versions will likely require you to change
how you use the project. Minor versions indicate new functionality that
doesn't impact the existing. Patch versions indicate
backwards-compatible fixes. Right now we are limiting this promise to
the command-line interface. A future release will introduce a stable and
mature Python API to the khmer project and at that time we will extend
the version system to include that API.

New items of note:
------------------

CITATION: Each script now outputs information on how to cite it. There
is a new paper to describes the project overall: MR Crusoe et al., 2014.
doi: 10.6084/m9.figshare.979190

The documentation for the scripts has undergone an overhaul. The scripts
now output extensive notes and the formal documentation website is
generated from the scripts themselves and will never be out of sync.

https://khmer.readthedocs.org/en/latest/scripts.html

Notable bugs fixed/issues closed:
---------------------------------

git clone of the khmer repo reqs > 0.5 GiB #223 @mr-c new khmer/file
module #357 @RamRS Floating point exception in count-overlap.py #282
@qingpeng add documentation for sample-reads-randomly #192 @mr-c only
build zlib and bzip2 when needed #168 @mr-c

Minor updates
-------------

khmer tools should output intelligent error messages when fed empty
files #135 @RamRS set IParser::ParserState::ParserState:fill\_id to zero
at initialization #356 @mr-c demote nose & sphinx to extra dependencies.
#351 @mr-c CID 1054792 (Medium) Uninitialized scalar field
(UNINIT\_CTOR) #179 @mr-c CID 1077117 (Medium): Division or modulo by
zero (DIVIDE\_BY\_ZERO) #182 @mr-c if --savehash is specified then don't
continue if there is not enough free disk space #245 @RamRS finish
fixing implicit downcasts #330 @mr-c Clean up compile warnings in
subset.cc #172 @mr-c all scripts need to output their version #236 @mr-c
environmental variables need documenting #303 @mr-c C++ code should be
consistently formatted #261 @mr-c Clean up ancillary files #146 @mr-c
squash option not implemented in abundance-dist-single.py #271 @RamRS
Add documentation on how to tie into a particular tagged version #29
@mr-c pip install -e fails with compile error #352 @mr-c remove the
unused KTable object #337 @luizirber zlib 1.2.3 -> zlib 1.2.8 #336 @mr-c
CID 1173035: Uninitialized scalar field (UNINIT\_CTOR) #311 @mr-c CID
1153101: Resource leak in object (CTOR\_DTOR\_LEAK) #309 @mr-c remove
khmer::read\_parsers::IParser::ParserState::thread\_id #323 @mr-c
several modifications about count-overlap.py script #324 @qingpeng fixed
runscript to handle SystemExit #332 @ctb CID 1063852: Uninitialized
scalar field (UNINIT\_CTOR) #313 @mr-c [infrastructure] update to new
Doxyfile format, make version number autoupdate #315 @mr-c Removed an
extraneous using namespace khmer; in kmer.hh, #276 @fishjord Minimum and
recommended python version #94 @mr-c KmerCount class appears to be
unused #302 @mr-c If loadhash is specified in e.g. normalize-by-median,
don't complain about default hashsize parameters #117 @RamRS

Known Issues
------------

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/dib-lab/khmer/issues/266

Some FASTQ files confuse our parser when running with more than one
thread. For example, while using load-into-counting.py. If you
experience this then add "--threads=1" to your command line. This issue
is being tracked in https://github.com/dib-lab/khmer/issues/249

If your k-mer table (hashfile) gets truncated, perhaps from a full
filesystem, then our tools currently will get stuck. This is being
tracked in https://github.com/dib-lab/khmer/issues/247 and
https://github.com/dib-lab/khmer/issues/96 and
https://github.com/dib-lab/khmer/issues/246

Paired-end reads from Casava 1.8 currently require renaming for use in
normalize-by-median and abund-filter when used in paired mode. The
integration of a fix for this is being tracked in
https://github.com/dib-lab/khmer/issues/23

annotate-partitions.py only outputs FASTA even if given a FASTQ file.
This issue is being tracked in
https://github.com/dib-lab/khmer/issues/46

A user reported that abundance-dist-single.py fails with small files and
many threads. This issue is being tracked in
https://github.com/dib-lab/khmer/issues/75

Contributors
------------

@camillescott, @mr-c, @ctb, @luizirber, @RamRS, @qingpeng

[0] http://semver.org/
