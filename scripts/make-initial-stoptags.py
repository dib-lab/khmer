#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Find an initial set of highly connected k-mers, to save on repartitioning time.

% python scripts/make-initial-stoptags.py <base>
"""

import argparse
import khmer
from khmer.file_api import check_file_status, check_space

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

    # @RamRS: This might need some more work
    infiles = [graphbase + '.ht', graphbase + '.tagset']
    if args.stoptags:
        infiles.append(args.stoptags)
    for _ in infiles:
        check_file_status(_)

    check_space(infiles)

    print 'loading htable %s.ht' % graphbase
    htable = khmer.load_hashbits(graphbase + '.ht')

    # do we want to load stop tags, and do they exist?
    if args.stoptags:
        print 'loading stoptags from', args.stoptags
        htable.load_stop_tags(args.stoptags)

    print 'loading tagset %s.tagset...' % graphbase
    htable.load_tagset(graphbase + '.tagset')

    ksize = htable.ksize()
    counting = khmer.new_counting_hash(ksize, args.min_hashsize, args.n_hashes)

    # divide up into SUBSET_SIZE fragments
    divvy = htable.divide_tags_into_subsets(args.subset_size)

    # pick off the first one
    if len(divvy) == 1:
        start, end = 0, 0
    else:
        start, end = divvy[:2]

    # partition!
    print 'doing pre-partitioning from', start, 'to', end
    subset = htable.do_subset_partition(start, end)

    # now, repartition...
    print 'repartitioning to find HCKs.'
    htable.repartition_largest_partition(subset, counting,
                                         EXCURSION_DISTANCE,
                                         EXCURSION_KMER_THRESHOLD,
                                         EXCURSION_KMER_COUNT_THRESHOLD)

    print 'saving stop tags'
    htable.save_stop_tags(graphbase + '.stoptags')

if __name__ == '__main__':
    main()
