#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Find an initial set of highly connected k-mers, to save on repartitioning time.

% python scripts/make-initial-stoptags.py <base>
"""

import sys
import argparse
import os
import khmer

#  Import fileapi from sandbox - temporary arrangement
current_file_path = os.path.realpath(__file__)
current_folder = os.path.dirname(current_file_path)
parent_folder = os.path.dirname(current_folder)
sandbox_folder = os.path.join(parent_folder, 'sandbox')
sys.path.append(sandbox_folder)

import fileApi

DEFAULT_SUBSET_SIZE = int(1e4)
DEFAULT_COUNTING_HT_SIZE = 3e6                # number of bytes
DEFAULT_COUNTING_HT_N = 4                     # number of counting hash tables

# Lump removal parameters.  Probably shouldn't be changed, but who knows?
#
# explanation:
#
# We will walk EXCURSION_DISTANCE out from each tag; if we find more than
# EXCURSION_KMER_THRESHOLD kmers within that range, this will be a "big"
# excursion and we will track all k-mers visited.  If we find that any
# k-mer has been visited more than EXCURSION_KMER_COUNT_THRESHOLD times,
# we will mark it as BAD and make it a stop tag for traversal.

# don't change these!
EXCURSION_DISTANCE = 40
EXCURSION_KMER_THRESHOLD = 200
EXCURSION_KMER_COUNT_THRESHOLD = 5

#


def main():
    parser = argparse.ArgumentParser(
        description="Find an initial set of highly connected k-mers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--n_hashes', '-N', type=int, dest='n_hashes',
                        default=DEFAULT_COUNTING_HT_N,
                        help='number of counting hash tables to use')
    parser.add_argument('--hashsize', '-x', type=float, dest='min_hashsize',
                        default=DEFAULT_COUNTING_HT_SIZE,
                        help='lower bound on counting hashsize to use')
    parser.add_argument('--subset-size', '-s', default=DEFAULT_SUBSET_SIZE,
                        dest='subset_size', type=float,
                        help='Set subset size (default 1e4 is prob ok)')
    parser.add_argument('--stoptags', '-S', dest='stoptags', default='',
                        help="Use stoptags in this file during partitioning")

    parser.add_argument('graphbase')

    args = parser.parse_args()

    graphbase = args.graphbase
    
    # Check input files exist
    infiles=[graphbase + '.ht', graphbase + '.tagset'] # @RamRS: This might need some more work
    for f in infiles:
        fileApi.check_file_status(f)

    # Check disk space availability
    freeSpace = fileApi.check_space(infiles)

    print 'loading ht %s.ht' % graphbase
    ht = khmer.load_hashbits(graphbase + '.ht')

    # do we want to load stop tags, and do they exist?
    if args.stoptags:
        print 'loading stoptags from', args.stoptags
        ht.load_stop_tags(args.stoptags)

    print 'loading tagset %s.tagset...' % graphbase
    ht.load_tagset(graphbase + '.tagset')

    K = ht.ksize()
    counting = khmer.new_counting_hash(K, args.min_hashsize, args.n_hashes)

    # divide up into SUBSET_SIZE fragments
    divvy = ht.divide_tags_into_subsets(args.subset_size)

    # pick off the first one
    if len(divvy) == 1:
        start, end = 0, 0
    else:
        start, end = divvy[:2]

    # partition!
    print 'doing pre-partitioning from', start, 'to', end
    subset = ht.do_subset_partition(start, end)

    # now, repartition...
    print 'repartitioning to find HCKs.'
    ht.repartition_largest_partition(subset, counting,
                                     EXCURSION_DISTANCE,
                                     EXCURSION_KMER_THRESHOLD,
                                     EXCURSION_KMER_COUNT_THRESHOLD)

    print 'saving stop tags'
    ht.save_stop_tags(graphbase + '.stoptags')

if __name__ == '__main__':
    main()
