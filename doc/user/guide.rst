.. vim: set filetype=rst

An assembly handbook for khmer - rough draft
############################################

:date: 2012-11-2

An increasing number of people are asking about using our assembly
approaches for things that we haven't yet written (or posted) papers
about.  Moreover, our assembly strategies themselves are also under
constant evolution as we do more research and find ever-wider
applicability of our approaches.

Note, this is an exact copy of `Titus' blog post, here
<http://ivory.idyll.org/blog/an-assembly-handbook-for-khmer.html>`__
-- go check the bottom of that for comments.

Authors
~~~~~~~

This handbook distills the cumulative expertise of Adina Howe, Titus
Brown, Erich Schwarz, Jason Pell, Camille Scott, Elijah Lowe, Kanchan
Pavangadkar, Likit Preeyanon, and others.

Introduction
~~~~~~~~~~~~

khmer is a general `framework for low-memory k-mer counting, filtering,
and advanced trickery <http://khmer.readthedocs.org/en/latest/>`__.

The latest source is always available `here
<https://github.com/dib-lab/khmer>`__.

khmer is really focused on short read data, and, more specifically,
Illumina, because that's where we have a too-much-data problem.
However, a lot of the prescriptions below can be adapted to longer
read technologies such as 454 and Ion Torrent without much effort.

Don't try to use our k-mer approaches with PacBio -- the error rate is
too high.

There are currently two papers available on khmer: the `partitioning
paper
<http://pnas.org/content/early/2012/07/25/1121464109.abstract>`__ and
the `digital normalization paper <http://arxiv.org/abs/1203.4802>`__.

There are many blog posts about this stuff on `Titus Brown's blog
<http://ivory.idyll.org/blog/>`__.  We will try to link them in where
appropriate.

Asking for help
~~~~~~~~~~~~~~~

There's some documentation here:

   https://khmer.readthedocs.org/en/latest/

There's also a khmer mailing list at lists.idyll.org that you can use to
get help with khmer.  To sign up, just go to 
`the khmer lists page <http://lists.idyll.org/listinfo/khmer>`__ and
subscribe.

Preparing your sequences
~~~~~~~~~~~~~~~~~~~~~~~~

Do all the quality filtering, trimming, etc. that you think you should do.

Most of the khmer tools currently work "out of the box" on interleaved
paired-end data.  Ask on the list if you're not sure.

All of our scripts will take in .fq or .fastq files as FASTQ, and all
other files as FASTA.  gzip files are always accepted.  Let us know if
not; that's a bug!

Most scripts *output* FASTA, and some mangle headers.  Sorry.  We're
working on outputting FASTQ for FASTQ input, and removing any header
mangling.

Picking k-mer table sizes and k parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For k-mer table sizes, read :doc:`choosing-table-sizes`

For k-mer sizes, we recommend k=20 for digital normalization and k=32
for partitioning; then assemble with a variety of k parameters.

Genome assembly, including MDA samples and highly polymorphic genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Apply digital normalization as follows.

Broadly, normalize each insert library separately, in the following way:

For high-coverage libraries (> ~50x), do three-pass digital
normalization: run normalize-by-median to C=20 and then run
filter-abund with C=1.  Now split out the remaining
paired-end/interleaved and single-end reads using
strip-and-split-for-assembly, and normalize-by-median the paired-end and
single-end files to C=5 (in that order).

For low-coverage libraries (< 50x) do single-pass digital normalization:
run normalize-by-median to C=10.

2. Extract any remaining paired-end reads and lump remaining orphan
reads into singletons using strip-and-split-for-assembly

3. Then assemble as normal, with appropriate insert size specs
etc. for the paired end reads.

You can read about this process in the `digital normalization paper
<http://arxiv.org/abs/1203.4802>`__.

mRNAseq assembly
~~~~~~~~~~~~~~~~

1. Apply single-pass digital normalization.

Run normalize-by-median to C=20.

2. Extract any remaining paired-end reads and lump remaining orphan
reads into singletons using strip-and-split-for-assembly

3. Then assemble as normal, with appropriate insert size specs
etc. for the paired end reads.

You can read about this process in the `digital normalization paper
<http://arxiv.org/abs/1203.4802>`__.

Metagenome assembly
~~~~~~~~~~~~~~~~~~~

1. Apply single-pass digital normalization.

Run normalize-by-median to C=20 (we've also found C=10 works fine).

2. Run filter-below-abund with C=50 (if you diginormed to C=10) or
C=100 (if you diginormed to C=20);

3. Partition reads with load-graph, etc. etc.

4. Assemble groups as normal, extracting paired-end reads and lumping
remaining orphan reads into singletons using
strip-and-split-for-assembly.

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

1. Apply single-pass digital normalization.

Run normalize-by-median to C=20.

2. Extract any remaining paired-end reads and lump remaining orphan
reads into singletons using strip-and-split-for-assembly

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
diginorm the Illumina data to C=20.

Quantifying mRNAseq or metagenomes assembled with digital normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For now, khmer only deals with assembly! So: assemble.  Then, go back
to your original, unnormalized reads, and map those to your assembly
with e.g. bowtie.  Then count as you normally would :).

Philosophy of digital normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The basic philosophy of digital normalization is "load your most
valuable reads first."  Diginorm gets rid of redundancy iteratively,
so you are more likely to retain the first reads fed in; this means
you should load in paired end reads, or longer reads, first.

Iterative and independent normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use :option:`--loadtable` and :option:`--savetable` to do iterative
normalizations on multiple files in multiple steps. For example, break ::

  normalize-by-median.py [ ... ] file1.fa file2.fa file3.fa

into multiple steps like so::

  normalize-by-median.py [ ... ] --savetable file1.ct file1.fa
  normalize-by-median.py [ ... ] --loadtable file1.ct --savetable file2.ct file2.fa
  normalize-by-median.py [ ... ] --loadtable file2.ct --savetable file3.ct file3.fa

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

Validating and comparing assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

More here soon :).

.. Check/validate assembly - look at high abundance kmers.
.. @@error trimming
