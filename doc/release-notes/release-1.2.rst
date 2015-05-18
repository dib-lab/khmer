khmer v1.2 release notes
========================

This is the v1.2 release of khmer: minor new features and bug fixes. The
start of this release cycle coincided with the Mozilla Science Lab
Global Sprint 2014. We honor and thank the 19 new contributors
(including four Michigan State University undergraduates) who
volunteered their time to contribute!

Docs at: https://khmer.readthedocs.org/en/v1.2/

New items of note:
------------------

@mr-c and @ctb are proud to announce khmer's code of conduct
http://khmer.readthedocs.org/en/v1.2/dev/CODE\_OF\_CONDUCT.html #664 All
scripts list which files have been created during their execution #477
@bocajnotnef All scripts now only output status messages to STDERR
instead of STDOUT #626 @b-wyss docs/ a fairly major re-organization and
brand new developer docs @ctb @mr-c load-into-counting.py:
``--summary-info``: machine readable summary in JSON or TSV format #649
@kdmurray91 scripts/extract-partitions.py: added documentation for .dist
columns #516 @chuckpr Makefile: a new target
``make install-dependencies`` is useful for developers #539 @mr-c
Sandbox scripts have been cleaned up, or removed (see the
sandbox/README.rst for details) #589 @ctb

Notable bugs fixed/issues closed:
---------------------------------

do-partition.py's excessive spawning of threads fixed. #637
@camillescott Fixed unique k-mer count reporting in load-graph,
load-into-counting, and normalize-by-median. #562 @mr-c Clarified and
test the requirement for a 64-bit operating system #529 @Echelon9
Removed some of the broken multi-threading options #511 @majoras-masque
Fix table.get("wrong\_length\_string") gives core dump #585 @Echelon9
filter-abund lists parameters that it doesn't use #524 @jstapleton
Reduction of memory required to run the test suite #542 @leogargu BibTeX
included in CITATIONS #541 @HLWiencko

Additional fixes/features
-------------------------

delete ScoringMatrix::assign as it is unused #502 @RodPic Root all of
our C++ exceptions to a common base exception #508 @iglpdc deleted
KhmerError #503 @drlabratory normalize-by-median reporting output after
main loop exits, in case it hadn't been triggered #586 @ctb Many issues
discovered by cppcheck cleaned up #506 @brtaylor92 Developers have a new
Makefile target to autofix formatting: ``make format`` #612 @brtaylor92
normalize-by-median.py test coverage increased #361 @SherineAwad Several
unused functions were removed #599 @brtaylor92 Developer docs now link
to the stdc++ docs as appropriate #629 @mr-c Added tests for
non-sequential access to input files #644 @bocajnotnef Removed
khmer/theading\_args.py #653 @bocajnotnef Improved test for maximum k
value #658 @pgarland ReadParser no longer crashes if n\_threads = 0 #86
@jiarong

Known issues:
-------------

All of these are pre-existing.

Multithreaded reading will drop reads. This major issue has been present
for several khmer releases and was only found via a much larger test
case that we had been previously using. Credit to @camillescott.
Workaround: disable threading. The next release will fix this and the
other FAST[AQ] parsing issues.
https://github.com/dib-lab/khmer/issues/681

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/dib-lab/khmer/issues/266

Some FASTQ files confuse our parser when running with more than one
thread. For example, while using load-into-counting.py. If you
experience this then add "--threads=1" to your command line. This issue
is being tracked in https://github.com/dib-lab/khmer/issues/249

If your k-mer table is truncated on write, an error may not be reported;
this is being tracked in https://github.com/dib-lab/khmer/issues/443.
However, khmer will now (correctly) fail when trying to read a truncated
file (See #333).

Paired-end reads from Casava 1.8 currently require renaming for use in
normalize-by-median and abund-filter when used in paired mode. The
integration of a fix for this is being tracked in
https://github.com/dib-lab/khmer/issues/23

Some scripts only output FASTA even if given a FASTQ file. This issue is
being tracked in https://github.com/dib-lab/khmer/issues/46

A user reported that abundance-dist-single.py fails with small files and
many threads. This issue is being tracked in
https://github.com/dib-lab/khmer/issues/75

Contributors
------------

@mr-c, @ctb, \*@bocajnotnef, \*@Echelon9, \*@jlippi, \*@kdmurray91,
@qingpeng, \*@leogargu, \*@jiarong, \*@brtaylor92, \*@iglpdc,
@camillescott, \*@HLWiencko, \*@cowguru2000, \*@drlabratory,
\*@jstapleton, \*@b-wyss, \*@jgluck, @fishjord, \*@SherineAwad,
\*@pgarland, \*@majoras-masque, @chuckpr, \*@RodPic, @luizirber,
\*@jrherr

``*`` Denotes new contributor
