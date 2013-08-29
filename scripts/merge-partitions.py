#! /usr/bin/env python
"""
Merge multiple pmap files into a single one.

% python scripts/merge-partitions.py <base>

Load <base>.subset.*.pmap and merge into a single pmap file.  Final
merged pmap file will be in <base>.pmap.merged.
"""

import sys
import argparse
import glob
import os

import khmer

DEFAULT_K = 32


def main():
    parser = argparse.ArgumentParser(description="Merge pmap files.")

    parser.add_argument('--ksize', '-k', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('--keep-subsets', dest='remove_subsets',
                        default=True, action='store_false',
                        help='Keep individual subsets (default: False)')
    parser.add_argument('graphbase')
    args = parser.parse_args()

    output_file = args.graphbase + '.pmap.merged'
    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print 'loading %d pmap files (first one: %s)' % (len(pmap_files),
                                                     pmap_files[0])

    K = args.ksize
    ht = khmer.new_hashbits(K, 1, 1)

    for pmap_file in pmap_files:
        print 'merging', pmap_file
        ht.merge_subset_from_disk(pmap_file)

    print 'saving merged to', output_file
    ht.save_partitionmap(output_file)

    if args.remove_subsets:
        print 'removing pmap files'
        for pmap_file in pmap_files:
            os.unlink(pmap_file)

if __name__ == '__main__':
    main()
