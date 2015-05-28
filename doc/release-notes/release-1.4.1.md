# khmer v1.4.1 release notes

This is the v1.4.1 release of khmer. Due to the upcoming Python 3 compatibility 
in both khmer and Screed we need to modify the dependency between khmer and the
Screed library to be only the existing version 0.8, and not some future
version.

If you have khmer 1.4 installed then there is no benefit to upgrading; this
point release is to keep `pip install khmer` still working when we release the
next version of Screed with Python 3 support. The next version of khmer, v2.0,
will also have Python 3 support.

Documentation is at https://khmer.readthedocs.org/en/v1.4.1/ (no changes from
v1.4)

## Known issues:

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more memory
than it was configured for. This is being investigated in
https://github.com/ged-lab/khmer/issues/266

Some scripts only output FASTA even if given a FASTQ file. This issue is being
tracked in https://github.com/ged-lab/khmer/issues/46
