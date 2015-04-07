.. vim: set filetype=rst

Known Issues
============

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/ged-lab/khmer/issues/266

If your k-mer table is truncated on write, an error may not be reported; this
is being tracked in https://github.com/ged-lab/khmer/issues/443.
However, khmer will now (correctly) fail when trying to read a truncated file
(See #333).

Paired-end reads from Casava 1.8 currently require renaming for use in
normalize-by-median and abund-filter when used in paired mode. The
integration of a fix for this is being tracked in
https://github.com/ged-lab/khmer/issues/23

Some scripts only output FASTA even if given a FASTQ file. This issue
is being tracked in https://github.com/ged-lab/khmer/issues/46
