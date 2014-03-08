#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import os
import gc
import glob

TRAVERSE_ON_UNPART = True
STOP_BIG_TRAVERSALS = True

already_part = sys.argv[1]
new_to_part = sys.argv[2]
basename = os.path.basename(new_to_part)
pmap_filename = sys.argv[3]

# if not os.path.exists(already_part):
#    print '%s doesn\'t exist! dying.' % already_part
#    sys.exit(0)

# create a fake-ish ht; K matters, but not hashtable size.
ht = khmer.load_hashbits(already_part + '.ht')
ht.load_tagset(already_part + '.tagset')
ht.merge_subset_from_disk(pmap_filename)

# find singletons
n_singletons = ht.find_unpart(
    new_to_part, TRAVERSE_ON_UNPART, STOP_BIG_TRAVERSALS)
print 'found:', n_singletons

print 'saving', basename + '.unpart'
n_partitions = ht.output_partitions(new_to_part, basename + '.unpart')
print 'saving', basename + '.pmap'
ht.save_partitionmap(basename + '.pmap')

###

(n_partitions, n_singletons) = ht.count_partitions()

print 'output partitions:', n_partitions
print 'pmap partitions:', n_partitions
print 'singletons:', n_singletons
