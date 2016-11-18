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
