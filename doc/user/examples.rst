.. vim: set filetype=rst

A few examples
==============

See the 'examples' subdirectory for complete examples.

STAMPS data set
---------------

The 'stamps' data set is a fake metagenome-like data set containing
two species, mixed at a 10:1 ratio.  The source genomes are
in 'data/stamps-genomes.fa'.  The reads file is in 'data/stamps-reads.fa.gz',
and consists of 100-base reads with a 1% error rate.

The example shows how to construct k-mer abundance histograms, as well
as the effect of digital normalization and partitioning on the k-mer
abundance distribution.

See `the script for running everything
<https://github.com/dib-lab/khmer/blob/master/examples/stamps/do.sh>`__
and `the IPython Notebook
<http://nbviewer.ipython.org/urls/raw.github.com/dib-lab/khmer/master/examples/stamps%2520k-mer%2520distributions.ipynb>`__.

For an overall discussion and some slides to explain what's going on,
visit `the Web site for a 2013 HMP metagenome assembly webinar that
Titus Brown gave <http://ged.msu.edu/angus/2013-hmp-assembly-webinar/exploring-stamps-data.html>`__.

