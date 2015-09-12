.. raw:: html

   <!--
      This file is part of khmer, https://github.com/dib-lab/khmer/, and is
      Copyright (C) 2014 Michigan State University
      It is licensed under the three-clause BSD license; see LICENSE.
      Contact: khmer-project@idyll.org
      
      Redistribution and use in source and binary forms, with or without
      modification, are permitted provided that the following conditions are
      met:
      
       * Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
      
       * Redistributions in binary form must reproduce the above
         copyright notice, this list of conditions and the following
         disclaimer in the documentation and/or other materials provided
         with the distribution.
      
       * Neither the name of the Michigan State University nor the names
         of its contributors may be used to endorse or promote products
         derived from this software without specific prior written
         permission.
      
      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
      HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      
      Contact: khmer-project@idyll.org
   -->

khmer v1.0.1 release notes
==========================

This is bugfix release. Note: the installation instructions have been
slightly simplified.

https://khmer.readthedocs.org/en/v1.0.1/

New items of note:
------------------

This release successfully installs and passes its unit tests on Debian
6.0 "Squeeze", Debian 7.0 "Wheezy", Fedora 19, OS X 7 "Lion", OS X 8
"Mountain Lion", Red Hat Enterprise Linux 6, Scientific Linux 6, Ubuntu
10.04 LTS, and Ubuntu 12.04 LTS. Thanks to the `UW-Madison Build and
Test Lab <https://www.batlab.org/>`__ for their `testing
infrastructure <http://submit-1.batlab.org/nmi/results/details?runID=247153>`__.

Notable bugs fixed/issues closed:
---------------------------------

fixed thread hanging issue #406 @ctb Explicit python2 invocation #404
@mr-c MANIFEST.in,setup.py: fix to correct zlib packaging #365 @mr-c
fixed check\_space\_for\_hashtable to use args.n\_tables #382 @ctb Bug
fix: make-initial-stoptags.py error on missing .ht input file, actual
input file is .pt #391 @mr-c

Minor updates
-------------

include calc-best-assembly.py in v1.0.1 #409 @ctb updated
normalize-by-median documentation for loadtable #378 @ctb updated
diginorm for new FP rate info; corrected spelling error #398 @ctb Add
spellcheck to code review checklist. #397 @ctb

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

@mr-c, @ctb, @luizirber, @RamRS, @ctSkennerton
