#! /bin/bash
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#

#
# download test data from:
#
#      https://s3.amazonaws.com/public.ged.msu.edu/khmer/iowa-corn-50m.fa.gz
#

KHMER_PATH=$1
export PYTHONPATH=$KHMER_PATH/python

SCRIPTPATH=$KHMER_PATH/scripts

# the next command will create a '50m.ht' and a '50m.tagset',
# representing the de Bruijn graph
${SCRIPTPATH}/load-graph.py -k 32 -N 4 -x 12e9 50m iowa-corn-50m.fa.gz 

# this will then partition that graph. should take a while.
# update threads to something higher if you have more cores.
# this creates a bunch of files, 50m.subset.*.pmap
${SCRIPTPATH}/partition-graph.py --threads 4 -s 1e5 50m

# now, merge the pmap files into one big pmap file, 50m.pmap.merged
${SCRIPTPATH}/merge-partitions.py 50m

# next, annotate the original sequences with their partition numbers.
# this will create iowa-corn-50m.fa.gz.part
${SCRIPTPATH}/annotate-partitions.py 50m iowa-corn-50m.fa.gz

# now, extract the partitions in groups into 'iowa-corn-50m.groupNNNN.fa'
${SCRIPTPATH}/extract-partitions.py iowa-corn-50m iowa-corn-50m.fa.gz.part

# at this point, you can assemble the group files individually.  Note,
# however, that the last one them is quite big?  this is because it's
# the lump! yay!

# if you want to break up the lump, go through the partitioning bit
# on the group file, but this time with a twist:
mv iowa-corn-50m.group0007.fa corn-50m.lump.fa

# create graph,
${SCRIPTPATH}/load-graph.py -x 8e9 lump corn-50m.lump.fa

# create an initial set of stoptags to help in knot-traversal; otherwise,
# partitioning and knot-traversal (which is systematic) is really expensive.
${SCRIPTPATH}/make-initial-stoptags.py lump

# now partition the graph, using the stoptags file
${SCRIPTPATH}/partition-graph.py --stoptags lump.stoptags lump

# use the partitioned subsets to find the k-mers that nucleate the lump
${SCRIPTPATH}/find-knots.py -x 2e8 -N 4 lump

# remove those k-mers from the fasta files
${SCRIPTPATH}/filter-stoptags.py *.stoptags corn-50m.lump.fa

# now, reload the filtered data set in and partition again.
${SCRIPTPATH}/load-graph.py -x 8e9 lumpfilt corn-50m.lump.fa.stopfilt
${SCRIPTPATH}/partition-graph.py -T 4 lumpfilt
${SCRIPTPATH}/merge-partitions.py lumpfilt
${SCRIPTPATH}/annotate-partitions.py lumpfilt corn-50m.lump.fa.stopfilt
${SCRIPTPATH}/extract-partitions.py corn-50m-lump corn-50m.lump.fa.stopfilt.part

# and voila, after all that, you should now have your de-knotted lump in
# corn-50m-lump.group*.fa.  The *.group????.fa files can now be
# assembled individually by your favorite assembler.
