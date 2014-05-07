#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import os
import gc
import glob

K = 32

filename = sys.argv[1]
basename = os.path.basename(filename)
subset_filenames = sys.argv[2:]

if not os.path.exists(filename):
    print '%s doesn\'t exist! dying.' % filename
    sys.exit(0)

print '---'
print 'outputting partitioned file to', basename + '.part'
print 'merged pmap will be in', basename + '.pmap.merged'
print '---'

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.new_hashbits(K, 1, 1)

# load & merge
for subset_file in subset_filenames:
    print '<-', subset_file
    ht.merge_subset_from_disk(subset_file)

# save merged partitionmap
if len(subset_filenames) > 1:
    ht.save_partitionmap(basename + '.pmap.merged')

# partition!
n_partitions = ht.output_partitions(filename, basename + '.part')
(n_partitions, n_singletons) = ht.count_partitions()
print 'output partitions:', n_partitions
print 'pmap partitions:', n_partitions
print 'singletons:', n_singletons
