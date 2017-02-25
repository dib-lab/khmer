..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2012-2015 Michigan State University
   Copyright (C) 2015 The Regents of the University of California.
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

==========================
Setting khmer memory usage
==========================

If you look at the documentation for the scripts (:doc:`scripts`) you'll
see a :option:`-M <load-into-counting.py -M>` parameter that sets the maximum
memory usage for any script that uses k-mer counting tables or k-mer graphs.
What is this?

khmer uses a special data structure that lets it store counting tables
and k-mer graphs in very low memory; the trick is that you must fix
the amount of memory khmer can use before running it. (See `Pell et
al., 2012 <http://www.ncbi.nlm.nih.gov/pubmed/22847406>`__ and `Zhang
et al., 2014 <http://www.ncbi.nlm.nih.gov/pubmed/25062443>`__ for the
details.)  This is what the :option:`-M <load-into-counting.py -M>` parameter
does.

If you set it too low, khmer will warn you to set it higher at the end.
See below for some good choices for various kinds of data.

**Note for khmer 1.x users:** As of khmer 2.0, the :option:`-M
<load-into-counting.py -M>` parameter sets the
:option:`-N <load-into-counting.py -N>`/:option:`--n_tables
<load-into-counting.py --n_tables>` and :option:`-x <load-into-counting.py -x>`
/:option:`--max-tablesize <load-into-counting.py --max-tablesize>` parameters
automatically. You can still set these parameters directly if you wish.

The really short version
========================

There is no way (except for experience, rules of thumb, and intuition) to
know what this parameter should be up front.  So, use the maximum
available memory::

  -M 16G

for a machine with 16 GB of free memory, for example. The supported suffixes
for setting memory usage are K, M, G, and T for kilobyte, megabyte, gigabyte,
and terabyte, respectively.

The short version
=================

This parameter specifies the maximum memory usage of the primary data
structure in khmer, which is basically N big hash tables of size x.
The **product** of the number of hash tables and the size of the hash
tables specifies the total amount of memory used, which is what the
:option:`-M <load-into-counting.py -M>` parameter sets.

These tables are used to track k-mers.  If they are too small, khmer
will fail in various ways (and will complain), but there is no harm
in making it too large. So, **the absolute safest thing to do is to
specify as much memory as is available**.  Most scripts will inform
you of the total memory usage, and (at the end) will complain if it's
too small.

Life is a bit more complicated than this, however, because some scripts --
:program:`load-into-counting.py` and :program:`load-graph.py` -- keep
ancillary information that will consume memory beyond this table data
structure.  So if you run out of memory, decrease the table size.

Also see the rules of thumb, below.

The long version
=====================

khmer's scripts, at their heart, represents k-mers in a very memory
efficient way by taking advantage of two data structures, `Bloom
filters <http://en.wikipedia.org/wiki/Bloom_filter>`__ and `Count-Min
Sketches <http://en.wikipedia.org/wiki/Count%E2%80%93min_sketch>`__, that are
both *probabilistic* and *constant memory*.  The "probabilistic" part
means that there are false positives: the less memory you use, the
more likely it is that khmer will think that k-mers are present when
they are not, in fact, present.

Digital normalization (:program:`normalize-by-median.py` and
:program:`filter-abund.py`) uses the Count-Min Sketch data structure.

Graph partitioning (:program:`load-graph.py` etc.) uses the Bloom filter data
structure.

The practical ramifications of this are pretty cool.  For example,
your digital normalization is guaranteed not to increase in memory
utilization, and graph partitioning is estimated to be 10-20x more
memory efficient than any other de Bruijn graph representation.  And
hash tables (which is what Bloom filters and Count-Min Sketches use)
are really fast and efficient.  Moreover, the optimal memory size for
these primary data structures is dependent on the number of k-mers,
but not explicitly on the size of k itself, which is very unusual.

In exchange for this memory efficiency, however, you gain a certain
type of parameter complexity.  Unlike your more typical k-mer package
(like the Velvet assembler, or Jellyfish or Meryl or Tallymer), you
are either guaranteed not to run out of memory (for digital
normalization) or much less likely to do so (for partitioning).

The biggest problem with khmer is that there is a minimum hash number
and size that you need to specify for a given number of k-mers, and
you cannot confidently predict what it is before actually loading in
the data.  This, by the way, is also true for de Bruijn graph
assemblers and all the other k-mer-based software -- the final memory
usage depends on the total number of k-mers, which in turn depends on
the true size of your underlying genomic variation (e.g. genome or
transcriptome size), the number of errors, and the k-mer size you
choose (the k parameter) `[ see Conway & Bromage, 2011 ]
<http://www.ncbi.nlm.nih.gov/pubmed?term=21245053>`__.  **The number
of reads or the size of your data set is only somewhat correlated with
the total number of k-mers.** Trimming protocols, sequencing depth,
and polymorphism rates are all important factors that affect k-mer
count.

The bad news is that we don't have good ways to estimate total k-mer
count a priori, although we can give you some rules of thumb, below.
In fact, counting the total number of distinct k-mers is a somewhat
annoying challenge.  Frankly, we recommend *just guessing* instead of
trying to be all scientific about it.

The good news is that you can never give khmer too much memory!  k-mer
counting and set membership simply gets more and more accurate as you
feed it more memory.  (Although there may be performance hits from
memory I/O, e.g.  `see the NUMA architecture
<http://en.wikipedia.org/wiki/Non-Uniform_Memory_Access>`__.)  The
other good news is that khmer can measure the false positive rate (FPR)
and detect dangerously low memory conditions.  For partitioning, we
actually *know* what a too-high FPR is -- our `k-mer
percolation paper <http://arxiv.org/abs/1112.4193>`__ lays out the
math.  For digital normalization, we assume that a FPR
of 20% is bad.  In both cases the data-loading scripts will exit with
an error-code.

If you insist on optimizing memory usage, the :program:`unique-kmers.py`
script will compute the approximate number of k-mers in a data set
fairly quickly. This number can be provided to several scripts via the
:option:`-U <load-into-counting.py -U>` option, which will use it to
calculate the FPR before processing any input data. If the amount of
requested memory yields an unacceptable FPR, the script will complain
loudly, giving you the chance to cancel the program before any time is
wasted. It will also report the minimum amount of memory required for
an acceptable FPR, so that you can immediately re-start the script with
the desired settings.

Rules of thumb
--------------

For digital normalization, we recommend:

 - ``-M 8G`` for any amount of sequencing for a single microbial genome,
   MDA-amplified or single colony.

 - ``-M 16G`` for up to a billion mRNAseq reads from any organism.  Past that,
   increase it.

 - ``-M 32G`` for most eukaryotic genome samples.

 - ``-M 32G`` will also handle most "simple" metagenomic samples (HMP on down)

 - For metagenomic samples that are more complex, such as soil or marine,
   start as high as possible.  For example, we are using ``-M 256G`` for
   ~300 Gbp of soil reads.

For partitioning of complex metagenome samples, we recommend starting
as high as you can -- something like half your system memory.  So if
you have 256 GB of RAM, use ``-M 128G`` which will use 128 GB of RAM
for the basic graph storage, leaving other memory for the ancillary
data structures.
