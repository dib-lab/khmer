khmer v1.4 release notes
========================

This is the v1.4 release of khmer featuring the results of our March and
April (PyCon) coding sprints and the 16 new contributors; the use of the
new v0.8 release of screed (the library we use for pure Python reading
of nucleotide sequence files); and the addition of @luizirber's
HyperLogLog counter for quick cardinality estimation.

Documentation is at https://khmer.readthedocs.org/en/v1.4/

New items of note:
------------------

Casava 1.8 read naming is now fully supported and in general the scripts
no longer mangle read names. Side benefits: ``split-paired-reads.py``
will no longer drop reads with 'bad' names; ``count-median.py`` can
generate output in CSV format. #759 #818 @ctb #873 @ahaerpfer

Most scripts now support a "broken" interleaved paired-read format for
FASTA/ FASTQ nucleotide sequence files.
```trim-low-abund.py`` <http://khmer.readthedocs.org/en/v1.4/user/scripts.html#trim-low-abund-py>`__
has been promoted from the sandbox as well (with streaming support).
#759 @ctb #963 @sguermond #933 @standage

The script to transform an interleaved paired-read nucleotide sequence
file into two files now allows one to name the output files which can be
useful in combination with named pipes for streaming processing #762
@ctb

Streaming everywhere: thanks to screed v0.8 we now support streaming of
almost all inputs and outputs. #830 @aditi9783 #812 @mr-c #917
@bocajnotnef #882 @standage

Need a quick way to count total number of unique k-mers in very low
memory? the ``unique-kmers.py`` script in the sandbox uses a HyperLogLog
counter to quickly (and with little memory) provide an estimate with a
controllable error rate. #257 #738 #895 #902 @luizirber

``normalize-by-median.py`` can now process both a paired interleaved
sequence file and a file of unpaired reads in the same invocation thus
removing the need to write the counting table to disk as required in the
workaround. #957 @susinmotion

Notable bugs fixed/issues closed:
---------------------------------

Paired-end reads from Casava 1.8 no longer require renaming for use in
``normalize-by-median.py`` and ``abund-filter.py`` when used in paired
mode #818 @ctb

Python version support clarified. We do not (yet) support Python 3.x
#741 @mr-c

If a single output file mode is chosen for normalize-by-median.py we now
default to overwriting the output. Appending the output is available by
using the append redirection operator from the shell. #843
@drtamermansour

Scripts that consume sequence data using C++ will now properly throw an
error on truncated files. #897 @kdmurray91 And while writing to disk we
properly check for errors #856 #962 @mr-c

``abundance-dist-single.py`` no longer fails with small files and many
threads. #900 @mr-c

Additional fixes/features
-------------------------

Of interest to users:
~~~~~~~~~~~~~~~~~~~~~

Many documentation updates #753 @PamelaM, #782 @bocajnotnef, #845
@alameldin, #804 @ctb, #870 @SchwarzEM, #953 #942 @safay,
#929,@davelin1, #687 #912 #926 @mr-c

Installation instructions for Conda, Arch Linux, and Mac Ports have been
added #723 @reedacartwright #952 @elmbeech #930 @ahaerpfer

The example script for the STAMPS database has been fixed to run
correctly #781 @drtamermansour

``split-paired-reads.py``: added ``-o`` option to allow specification of
an output directory #752 @bede

Fixed a string formatting and a boundry error in
``sample-reads-randomly.py`` #773 @qingpeng #995 @ctb

CSV output added to ``abundance-dist.py``, ``abundance-dist-single.py``,
and ``count-overlap.py``, and ``readstats.py`` #831 #854 #855
@drtamermansour #959 @anotherthomas

TSV/JSON output of ``load-into-counting.py`` enhanced with the total
number of reads processed #996 @kdmurray91 Output files are now also
checked to be writable *before* loading the input files #672 @pgarland
@bocajnotnef

``interleave-reads.py`` now prints the output filename nicely #827
@kdmurray91

Cleaned up error for input file not existing #772 @jessicamizzi #851
@ctb

Fixed error in ``find-knots.py`` #860 @TheOneHyer

The help text for ``load-into-counting.py`` for the
``--no-bigcounts``/``-b`` flag has been clarified #857 @kdmurray91

@lexnederbragt confirmed an old bug has been fixed with his test for
whitespace in sequence identifiers interacting with the
``extract-partitions.py`` script #979

Now safe to copy-and-paste from the user documentation as the smart
quotes have been turned off. #967 @ahaerpfer

The script ``make-coverage.py`` has been restored to the sandbox. #920
@SherineAwad

``normalize-by-median.py`` will warn if two of the input files have the
same name #932 @elmbeech

Of interest to developers:
~~~~~~~~~~~~~~~~~~~~~~~~~~

Switched away from using ``--user`` install for developers #740 @mr-c
@drtamermansour & #883 @standage

Developers can now see a summary of important Makefile targets via
``make help`` #783 @standage

The unused ``khmer.load_pe`` module has been removed #828 @kdmurray91

Versioneer bug due to new screed release was squashed #835 @mr-c

A Python 2.6 and 2.7.2 specific bug was worked around #869 @kdmurray91
@ctb

Added functions hash\_find\_all\_tags\_list and
hash\_get\_tags\_and\_positions to CountingHash objects #749 #765 @ctb

The ``make diff-cover`` and ChangeLog formatting requirements have been
added to checklist #766 @mr-c

A useful message is now presented if large tables fail to allocate
enough memory #704 @mr-c

A checklist for developers adding new CPython types was added #727 @mr-c

The sandbox graduation checklist has been updated to include streaming
support #951 @sguermond

Specific policies for sandbox/ and scripts/ content, and a process for
adding new command line scripts into scripts/ have been added to the
developer documentation #799 @ctb

Sandbox scripts update: corrected #! Python invocation #815 @Echelon9,
executable bits, copyright headers, no underscores in filenames #823
#826 #850 @alameldin several scripts deleted, docs + requirements
updated #852 @ctb

Avoid running big-memory tests on OS X #819 @ctb

Unused callback code was removed #698 @mr-c

The CPython code was updated to use the new checklist and follow
additional best practices #785 #842 @luizirber

Added a read-only view of the raw counting tables #671 @camillescott
#869 @kdmurray91

Added a Python method for quickly getting the number of underlying
tables in a counting or presence table #879 #880 @kdmurray91

The C++ library can now be built separately for the brave and curious
developer #788 @kdmurray91

The ReadParser object now keeps track of the number of reads processed
#877 @kdmurray91

Documentation is now reproducible #886 @mr-c

Python future proofing: specify floor division #863 @mr-c

Miscellaneous spelling fixes; thanks codespell! #867 @mr-c

Debian package list update #984 @mr-c

``khmer.kfile.check_file_status()`` has been renamed to
``check_input_files()`` #941 @proteasome ``filter-abund.py`` now uses it
to check the input counting table #931 @safay

``normalize-by-median.py`` was refactored to not pass the ArgParse
object around #965 @susinmotion

Developer communication has been clarified #969 @sguermond

Tests using the 'fail\_okay=true' parameter to ``runscript`` have been
updated to confirm the correct error occurred. 3 faulty tests were fixed
and the docs were clarified #968 #971 @susinmotion

FASTA test added for ``extract-long-sequences.py`` #901 @jessicamizzi

'added silly test for empty file warning' #557 @wltrimbl @bocajnotnef

A couple tests were made more resilient and some extra error checking
added in CPython land #889 @mr-c

Copyright added to pull request checklist #940 @sguermond

``khmer_exception``\ s are now based on ``std::string``\ s which plugs a
memory leak #938 @anotherthomas

Python docstrings were made PEP257 compliant #936 @ahaerpfer

Some C++ comments were converted to be Doxygen compliant #950
@josiahseaman

The counting and presence table warning logic was refactored and
centralized #944 @susinmotion

The release checklist was updated to better run the post-install tests
#911 @mr-c

The unused method ``find_all_tags_truncate_on_abundance`` was removed
from the CPython API #924 @anotherthomas

OS X warnings quieted #887 @mr-c

Known issues:
-------------

All of these are pre-existing.

Some users have reported that normalize-by-median.py will utilize more
memory than it was configured for. This is being investigated in
https://github.com/dib-lab/khmer/issues/266

Some scripts only output FASTA even if given a FASTQ file. This issue is
being tracked in https://github.com/dib-lab/khmer/issues/46

Contributors
------------

@ctb, @kdmurray91, @mr-c, @drtamermansour, @luizirber, @standage,
@bocajnotnef, \*@susinmotion, @jessicamizzi, \*@elmbeech,
\*@anotherthomas, \*@sguermond, \*@ahaerpfer, \*@alameldin,
\*@TheOneHyer, \*@aditi9783, \*@proteasome, \*@bede, \*@davelin1,
@Echelon9, \*@reedacartwright, @qingpeng, \*@SchwarzEM, \*@scottsievert,
@PamelaM, @SherineAwad, \*@josiahseaman, \*@lexnederbragt,

\* Indicates new contributors

Issue reporters
---------------

@moorepants, @teshomem, @macmanes, @lexnederbragt, @r-gaia-cs,
@magentashades
