..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2011-2015 Michigan State University
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

*********************
Introduction to khmer
*********************

Introduction
============

khmer is a software library and toolkit for k-mer based analysis and transformation of nucleotide sequence data.
The primary focus of development has been scaling assembly of metagenomes and transcriptomes.

khmer supports a number of transformations, both inexact transformations (abundance filtering; error trimming) and exact transformations (graph-size filtering, to throw away disconnected reads; partitioning, to split reads into disjoint sets).
All of these transformations operate with constant memory consumption (with the exception of partitioning), and typically require less memory than is required to assemble the data.

Most of khmer is built around a single underlying probabilistic data structure known as a `Bloom filter <http://en.wikipedia.org/wiki/Bloom_filter>`__ (also see `Count-Min Sketch <http://dimacs.rutgers.edu/~graham/pubs/papers/cm-full.pdf>`__ and `These Are Not The k-mers You're Looking For <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4111482/>`__).
In khmer this data structure is implemented as a set of hash tables, each of different size, with no collision detection.
These hash tables can store either the presence (Bloom filter) or frequency (Count-Min Sketch) of specific k-mers.
The lack of collision detection means that the data structure may report a k-mer as being *present* when it is not, in fact, in the data set.
However, it will never incorrectly report a k-mer as being absent when it *truly is* present.
This one-sided error makes the Bloom filter very useful for certain kinds of operations.

khmer supports arbitrarily large k-sizes, although certain graph-based operations are limited to k <= 32.

The khmer core library is implemented in C++, while all of the khmer scripts and tests access the core library via a Python wrapper.

Tutorials highlighting khmer are available at `khmer-protocols <http://khmer-protocols.readthedocs.io>`__ and `khmer-recipes <http://khmer-recipes.readthedocs.io>`__.
The former provides detailed protocols for using khmer to analyze either a transcriptome or a metagenome.
The latter provides individual recipes for using khmer in a variety of sequence-oriented tasks such as extracting reads by coverage, estimating a genome or metagenome size from unassembled reads, and error-trimming reads via streaming k-mer abundance.

Using khmer
===========

khmer comes "out of the box" with a number of scripts that make it immediately useful for a few different operations, including (but not limited to) the following.

 - normalizing read coverage ("digital normalization")
 - dividing reads into disjoint sets that do not connect ("partitioning")
 - eliminating reads that will not be used by a de Bruijn graph assembler;
 - removing reads with low- or high-abundance k-mers;
 - trimming reads of certain kinds of sequencing errors;
 - counting k-mers and estimating data set coverage based on k-mer counts;
 - running Velvet and calculating assembly statistics;
 - converting FASTQ to FASTA;
 - converting between paired and interleaved formats for paired FASTQ data

Practical considerations
========================

The most important thing to think about when using khmer is whether or not the transformation or filter you're applying is appropriate for the data you're trying to assemble.
For example, two of the most powerful operations available in khmer, graph-size filtering and graph partitioning, only make sense for assembly datasets with many theoretically unconnected components.
This is typical of metagenomic data sets.
Also, while digital normalization can be helpful for transcriptome *assembly*, it is inappropriate for other RNA-seq applications, such as differential expression analysis, that rely on signal from variable coverage.

The second most important consideration is memory usage.
The effectiveness of all of the Bloom filter-based functions (which is everything interesting in khmer!) depends critically on having enough memory to do a good job.
See :doc:`user/choosing-table-sizes` for more information.

Copyright and license
=====================

Portions of khmer are Copyright California Institute of Technology, where the exact counting code was first developed.
All other code developed through 2014 is copyright Michigan State University.
Portions are copyright Michigan State University and Regents of the University of California.
All the code is freely available for use and re-use under the BSD License.
Distribution, modification and redistribution, incorporation into other software, and pretty much everything else is allowed.
