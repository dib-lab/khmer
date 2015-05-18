khmer v1.3 release notes
========================

This is the v1.3 release of khmer featuring a new FAST[AQ] parser from
the SeqAn project.

Docs at: https://khmer.readthedocs.org/en/v1.3/

New items of note:
------------------

Fixes the two multithreaded reading of sequence files issues: FASTQ
parsing and the recently found read dropping issue. Several khmer
scripts now support reading from non-seekable plain and gziped FAST[AQ]
files (a.k.a pipe or streaming support). @mr-c #642

Notable bugs fixed/issues closed:
---------------------------------

restore threading to load-graph.py #699 @mr-c

Additional fixes/features
-------------------------

increase filter\_abund.py coverage #568 @wrightmhw Provide scripts/
testing coverage for check\_space\_for\_hashtable #386 #678 #718 @b-wyss
Use absolute URI in CODE\_OF\_CONDUCT #684 @jsspencer give SeqAn credit
#712 @mr-c Added testing to make sure all sandbox scripts are
import-able and execfile-able. #709 @ctb reduce memory requirements to
run tests #701 @ctb Two minor bug fixes to sandbox scripts #706 @ctb
Upgrade of trim-low-abund for better, more profitable streaming. #601
@ctb Add --force or --expert or --ignore flag to all khmer scripts that
do sanity checking #399 #647 @jessicamizzi Add XDECREF for returned read
tuple in ReadParser.read\_pair\_iterator() #693 @mr-c @camillescott

Known issues:
-------------

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in #266

If your k-mer table is truncated on write, an error may not be reported;
this is being tracked in https://github.com/dib-lab/khmer/issues/443.
However, khmer will now (correctly) fail when trying to read a truncated
file (See #333).

Paired-end reads from Casava 1.8 currently require renaming for use in
normalize-by-median and abund-filter when used in paired mode. The
integration of a fix for this is being tracked in #23

Some scripts only output FASTA even if given a FASTQ file. This issue is
being tracked in #46

A user reported that abundance-dist-single.py fails with small files and
many threads. This issue is being tracked in #75

Contributors
------------

@mr-c, @ctb, @camillescott, @b-wyss, @wrightmhw, @jsspencer
