#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Find highly-connected k-mers.

k-mers are output into a .stoptags file, for later use in partitioning.

% python scripts/find-knots.py <base>
"""
from __future__ import print_function

import argparse
import glob
import os
import textwrap
import khmer
import sys
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info

# counting hash parameters.
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
EXCURSION_KMER_COUNT_THRESHOLD = 2
# EXCURSION_KMER_COUNT_THRESHOLD=5 # -- works ok for non-diginormed data


def get_parser():
    epilog = """
    Load an k-mer presence table/tagset pair created by load-graph, and a set
    of pmap files created by partition-graph. Go through each pmap file,
    select the largest partition in each, and do the same kind of traversal as
    in :program:`make-initial-stoptags.py` from each of the waypoints in that
    partition; this should identify all of the HCKs in that partition. These
    HCKs are output to <graphbase>.stoptags after each pmap file.

    Parameter choice is reasonably important. See the pipeline in
    :doc:`partitioning-big-data` for an example run.

    This script is not very scalable and may blow up memory and die horribly.
    You should be able to use the intermediate stoptags to restart the
    process, and if you eliminate the already-processed pmap files, you can
    continue where you left off.
    """
    parser = argparse.ArgumentParser(
        description="Find all highly connected k-mers.",
        epilog=textwrap.dedent(epilog),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--n_tables', '-N', type=int,
                        default=DEFAULT_COUNTING_HT_N,
                        help='number of k-mer counting tables to use')
    parser.add_argument('--min-tablesize', '-x', type=float,
                        default=DEFAULT_COUNTING_HT_SIZE, help='lower bound on'
                        ' the size of the k-mer counting table(s)')
    parser.add_argument('graphbase', help='Basename for the input and output '
                        'files.')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Continue past warnings')
    return parser


def main():
    info('find-knots.py', ['graph'])
    args = get_parser().parse_args()

    graphbase = args.graphbase

    # @RamRS: This might need some more work
    infiles = [graphbase + '.pt', graphbase + '.tagset']
    if os.path.exists(graphbase + '.stoptags'):
        infiles.append(graphbase + '.stoptags')
    for _ in infiles:
        check_input_files(_, args.force)

    check_space(infiles, args.force)

    print('loading k-mer presence table %s.pt' % graphbase, file=sys.stderr)
    htable = khmer.load_hashbits(graphbase + '.pt')

    print('loading tagset %s.tagset...' % graphbase, file=sys.stderr)
    htable.load_tagset(graphbase + '.tagset')

    initial_stoptags = False    # @CTB regularize with make-initial
    if os.path.exists(graphbase + '.stoptags'):
        print('loading stoptags %s.stoptags' % graphbase, file=sys.stderr)
        htable.load_stop_tags(graphbase + '.stoptags')
        initial_stoptags = True

    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print('loading %d pmap files (first one: %s)' %
          (len(pmap_files), pmap_files[0]), file=sys.stderr)
    print('---', file=sys.stderr)
    print('output stoptags will be in',
          graphbase + '.stoptags', file=sys.stderr)
    if initial_stoptags:
        print(
            '(these output stoptags will include the already-loaded set)',
            file=sys.stderr)
    print('---', file=sys.stderr)

    # create counting hash
    ksize = htable.ksize()
    counting = khmer.new_counting_hash(ksize, args.min_tablesize,
                                       args.n_tables)

    # load & merge
    for index, subset_file in enumerate(pmap_files):
        print('<-', subset_file, file=sys.stderr)
        subset = htable.load_subset_partitionmap(subset_file)

        print('** repartitioning subset... %s' % subset_file, file=sys.stderr)
        htable.repartition_largest_partition(subset, counting,
                                             EXCURSION_DISTANCE,
                                             EXCURSION_KMER_THRESHOLD,
                                             EXCURSION_KMER_COUNT_THRESHOLD)

        print('** merging subset... %s' % subset_file, file=sys.stderr)
        htable.merge_subset(subset)

        print('** repartitioning, round 2... %s' %
              subset_file, file=sys.stderr)
        size = htable.repartition_largest_partition(
            None, counting, EXCURSION_DISTANCE, EXCURSION_KMER_THRESHOLD,
            EXCURSION_KMER_COUNT_THRESHOLD)

        print('** repartitioned size:', size, file=sys.stderr)

        print('saving stoptags binary', file=sys.stderr)
        htable.save_stop_tags(graphbase + '.stoptags')
        os.rename(subset_file, subset_file + '.processed')
        print('(%d of %d)\n' % (index, len(pmap_files)), file=sys.stderr)

    print('done!', file=sys.stderr)

if __name__ == '__main__':
    main()
