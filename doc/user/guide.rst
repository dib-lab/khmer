..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2011-2015 Michigan State University
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

An assembly handbook for khmer - rough draft
############################################

:date: 2012-11-2

An increasing number of people are asking about using our assembly
approaches for things that we haven't yet written (or posted) papers
about.  Moreover, our assembly strategies themselves are also under
constant evolution as we do more research and find ever-wider
applicability of our approaches.

Note, this is modified from `Titus' blog post, here
<http://ivory.idyll.org/blog/an-assembly-handbook-for-khmer.html>`__
-- go check the bottom of that for comments.

Authors
~~~~~~~

This handbook distills the cumulative expertise of Adina Howe, Titus
Brown, Erich Schwarz, Jason Pell, Camille Scott, Elijah Lowe, Kanchan
Pavangadkar, Likit Preeyanon, and others.

Introduction
~~~~~~~~~~~~

khmer is really focused on short read data, and, more specifically,
Illumina, because that's where we have a too-much-data problem.
However, a lot of the prescriptions below can be adapted to longer
read technologies such as 454 and Ion Torrent without much effort.

Don't try to use our k-mer approaches with PacBio -- the error rate is
too high.

There are many blog posts about this stuff on `Titus Brown's blog
<http://ivory.idyll.org/blog/>`__.  We will try to link them in where
appropriate.

Preparing your sequences
~~~~~~~~~~~~~~~~~~~~~~~~

Do all the quality filtering, trimming, etc. that you think you should do.

The khmer tools work "out of the box" on interleaved paired-end data.

All of our scripts will take in .fq or .fastq files as FASTQ, and all
other files as FASTA.  gzip files are always accepted.  Let us know if
not; that's a bug!

Genome assembly, including MDA samples and highly polymorphic genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Apply digital normalization as follows.

Broadly, normalize each insert library separately, in the following way:

For high-coverage libraries (> ~50x), do three-pass digital
normalization: run :program:`normalize-by-median.py` with :option:`--cutoff=20
<normalize-by-median.py --cutoff>` and then run :program:`filter-abund.py` with
:option:`--cutoff=2 <filter-abund.py --cutoff>`.  Now split out the remaining
paired-end/interleaved and single-end reads using
:program:`extract-paired-reads.py`, and run :program:`normalize-by-median.py`
on the paired-end and single-end files (using :option:`--unpaired-reads
<normalize-by-median.py --unpaired-reads>`) with :option:`--cutoff=5
<normalize-by-median.py --cutoff>`.

For low-coverage libraries (< 50x) do single-pass digital normalization:
run :program:`normalize-by-median.py` to :option:`--cutoff=10
<normalize-by-median.py --cutoff>`.

2. Extract any remaining paired-end reads and lump remaining orphan
   reads into singletons using :program:`extract-paired-reads.py`

3. Then assemble as normal, with appropriate insert size specs
   etc. for the paired end reads.

You can read about this process in the `digital normalization paper
<http://arxiv.org/abs/1203.4802>`__.

mRNAseq assembly
~~~~~~~~~~~~~~~~

1. Apply single-pass digital normalization.
   Run :program:`normalize-by-median.py` with :option:`--cutoff=20
   <normalize-by-median.py --cutoff>`.

2. Extract any remaining paired-end reads and lump remaining orphan
   reads into singletons using :program:`extract-paired-reads.py`

3. Then assemble as normal, with appropriate insert size specs
   etc. for the paired end reads.

You can read about this process in the `digital normalization paper
<http://arxiv.org/abs/1203.4802>`__.

Metagenome assembly
~~~~~~~~~~~~~~~~~~~

1. Apply single-pass digital normalization.
   Run :program:`normalize-by-median.py` with :option:`--cutoff=20
   <normalize-by-median.py --cutoff>` (we've also found :option:`--cutoff=10
   <normalize-by-median.py --cutoff>` works
   fine).

2. Run ``sandbox/filter-below-abund.py`` with ``--cutoff=50`` (if you
   ran :program:`normalize-by-median.py` with :option:`--cutoff=10
   <normalize-by-median.py --cutoff>`) or wiht ``--cutoff=100`` if you ran
   :program:`normalize-by-median.py` with :option:`--cutoff=20
   <normalize-by-median.py --cutoff>`)

3. Partition reads with :program:`load-graph.py`, etc. etc.

4. Assemble groups as normal, extracting paired-end reads and lumping
   remaining orphan reads into singletons using
   :program:`extract-paired-reads.py`.

(We actually use Velvet at this point, but there should be no harm in
using a metagenome assembler such as MetaVelvet or MetaIDBA or
SOAPdenovo.)

Read more about this in the `partitioning
<http://pnas.org/content/early/2012/07/25/1121464109.abstract>`__
paper.  We have some upcoming papers on partitioning and metagenome
assembly, too; we'll link those in when we can.

Metatranscriptome assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~

(Not tested by us!)

1. Apply single-pass digital normalization by running
   :program:`normalize-by-median.py` with :option:`--cutoff=20
   <normalize-by-median.py --cutoff>`.

2. Extract any remaining paired-end reads and lump remaining orphan
   reads into singletons using :program:`extract-paired-reads.py`.

3. Then assemble with a genome or metagenome assembler, *not* an
   mRNAseq assembler. Use appropriate insert size specs etc. for the
   paired end reads.

Preprocessing Illumina for other applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(Not tested by us!)

Others have told us that you can apply digital normalization to
Illumina data prior to using Illumina for `RNA scaffolding
<http://www.ncbi.nlm.nih.gov/pubmed?term=20980554>`__ or `error
correcting PacBio reads
<http://www.ncbi.nlm.nih.gov/pubmed?term=22750884>`__.

Our suggestion for this, based on no evidence whatsoever, is to
run :program:`normalize-by-median.py` with :option:`--cutoff=20
<normalize-by-median.py --cutoff>` on the Illumina data.

Quantifying mRNAseq or metagenomes assembled with digital normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For now, khmer only deals with assembly! So: assemble.  Then, go back
to your original, unnormalized reads, and map those to your assembly
with e.g. bowtie.  Then count as you normally would).

Philosophy of digital normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The basic philosophy of digital normalization is "load your most
valuable reads first."  Diginorm gets rid of redundancy iteratively,
so you are more likely to retain the first reads fed in; this means
you should load in paired end reads, or longer reads, first.

Iterative and independent normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use :option:`--loadgraph <normalize-by-median.py --loadgraph>` and
:option:`--savegraph <normalize-by-median.py --savegraph>` to do iterative
normalizations on multiple files in multiple steps. For example, break ::

  normalize-by-median.py [ ... ] file1.fa file2.fa file3.fa

into multiple steps like so::

  normalize-by-median.py [ ... ] --savegraph file1.ct file1.fa
  normalize-by-median.py [ ... ] --loadgraph file1.ct --savegraph file2.ct file2.fa
  normalize-by-median.py [ ... ] --loadgraph file2.ct --savegraph file3.ct file3.fa

The results should be identical!

If you want to independently normalize multiple files for speed reasons, go
ahead.  Just remember to do a combined normalization at the end.  For example,
instead of ::

  normalize-by-median.py [ ... ] file1.fa file2.fa file3.fa

you could do ::

  normalize-by-median.py [ ... ] file1.fa
  normalize-by-median.py [ ... ] file2.fa
  normalize-by-median.py [ ... ] file3.fa

and then do a final ::

  normalize-by-median.py [ ... ] file1.fa.keep file2.fa.keep file3.fa.keep

The results will not be identical, but should not differ
significantly.  The multipass approach will take more total time but
may end up being faster walltime because you can execute the
independent normalizations on multiple computers.

For a cleverer approach that we will someday implement, read `the
Beachcomber's Dilemma
<http://ivory.idyll.org/blog/beachcombers-dilemma.html>`__.

.. Validating and comparing assemblies
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. More here soon :).

.. Check/validate assembly - look at high abundance kmers.
.. @@error trimming
