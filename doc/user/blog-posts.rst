..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2010-2015 Michigan State University
   Copyright (C) 2015-2016 The Regents of the University of California.
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

=======================================
Blog posts and additional documentation
=======================================

Hashtable and filtering
=======================

The basic inexact-matching approach used by the hashtable code is
described in this blog post:

   http://ivory.idyll.org/blog/jul-10/kmer-filtering

A test data set (soil metagenomics, 88m reads, 10gb) is here:

   http://ci.oxli.org/userContent/88m-reads.fa.gz

Illumina read abundance profiles
================================

khmer can be used to look at systematic variations in k-mer statistics
across Illumina reads; see, for example, this blog post:

   http://ivory.idyll.org/blog/jul-10/illumina-read-phenomenology

The `fasta-to-abundance-hist
<http://github.com/ctb/khmer/blob/master/sandbox/fasta-to-abundance-hist.py>`__
and `abundance-hist-by-position
<http://github.com/ctb/khmer/blob/master/sandbox/abundance-hist-by-position.py>`__
scripts can be used to generate the k-mer abundance profile data, after
loading all the k-mer counts into a .ct file::

   # first, load all the k-mer counts:
   load-into-counting.py -k 20 -x 1e7 25k.ct data/25k.fq.gz

   # then, build the '.freq' file that contains all of the counts by position
   python sandbox/fasta-to-abundance-hist.py 25k.ct data/25k.fq.gz

   # sum across positions.
   python sandbox/abundance-hist-by-position.py data/25k.fq.gz.freq > out.dist

The hashtable method 'dump_kmers_by_abundance' can be used to dump
high abundance k-mers, but we don't have a script handy to do that yet.

You can assess high/low abundance k-mer distributions with the
`hi-lo-abundance-by-position script <http://github.com/ctb/khmer/blob/master/sandbox/hi-lo-abundance-by-position.py>`__::

   load-into-counting.py -k 20 25k.ct data/25k.fq.gz
   python sandbox/hi-lo-abundance-by-position.py 25k.ct data/25k.fq.gz

This will produce two output files, <filename>.pos.abund=1 and
<filename>.pos.abund=255.


Finding valleys/minima in *k*-mer abundance profiles
====================================================

Using *k*-mer abundance profiles to dynamically calculate the abundance threshold separating erroneous *k*-mers from real *k*-mers is described in this blog post:

    https://bitsandbugs.org/2016/07/29/mash-and-khmer-abundance/
