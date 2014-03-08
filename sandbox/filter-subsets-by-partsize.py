#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import khmer
import gc
import os

MIN_PARTITION_SIZE = 5          # 3 is the smallest reasonable threshold
K = 32


def main(subset_filenames):
    print 'K', K
    print 'MIN SIZE', MIN_PARTITION_SIZE
    print '--'

    # create an empty hashtable & load in the tags
    ht = khmer.new_hashbits(32, 1, 1)
    tagmap = ht.new_tagmap()

    # find the maximum partition size for each tag, across all subsets
    for filename in subset_filenames:
        print 'maxifying:', filename
        subset = ht.load_subset_partitionmap(filename)
        ht.subset_maxify_partition_size(subset, tagmap)
        del subset
        gc.collect()

    # filter tags based on the max partition size to which they belong
    print 'discarding'
    ht.discard_tags(tagmap, MIN_PARTITION_SIZE)

    # finally, filter each subset filename and save.
    for filename in subset_filenames:
        print 'loading x 2', filename
        subset = ht.load_subset_partitionmap(filename)
        print 'filtering', filename
        ht.subset_filter_against_tags(subset, tagmap)

        dir = os.path.dirname(filename)
        new_filename = 'filtered_' + os.path.basename(filename)
        new_filename = os.path.join(dir, new_filename)

        print 'saving', new_filename
        ht.save_subset_partitionmap(subset, new_filename)

        del subset
        gc.collect()

if __name__ == '__main__':
    subset_filenames = sys.argv[1:]
    main(subset_filenames)
