#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=invalid-name,missing-docstring
"""
Find highly-connected k-mers.

k-mers are output into a .stoptags file, for later use in partitioning.

% python scripts/find-knots.py <base>
"""

import glob
import os
import textwrap
import khmer
import sys
from khmer import Nodegraph, SubsetPartition
from khmer.kfile import check_input_files, check_space
from khmer import khmer_args
from khmer.khmer_args import (build_counting_args, sanitize_help)

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
    epilog = """\
    Load an k-mer nodegraph/tagset pair created by
    :program:`load-graph.py`, and a set of pmap files created by
    :program:`partition-graph.py`. Go through each pmap file,
    select the largest partition in each, and do the same kind of traversal as
    in :program:`make-initial-stoptags.py` from each of the waypoints in that
    partition; this should identify all of the Highly Connected Kmers in that
    partition. These HCKs are output to ``<graphbase>.stoptags`` after each
    pmap file.

    Parameter choice is reasonably important. See the pipeline in
    :doc:`partitioning-big-data` for an example run.

    This script is not very scalable and may blow up memory and die horribly.
    You should be able to use the intermediate stoptags to restart the
    process, and if you eliminate the already-processed pmap files, you can
    continue where you left off.
    """
    parser = build_counting_args(
        descr="Find all highly connected k-mers.",
        epilog=textwrap.dedent(epilog),
        citations=['graph'])

    parser.add_argument('graphbase', help='Basename for the input and output '
                        'files.')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Continue past warnings')
    return parser


def main():
    parser = get_parser()
    parser.epilog = parser.epilog.replace(
        ":doc:`partitioning-big-data`",
        "http://khmer.readthedocs.io/en/stable/user/"
        "partitioning-big-data.html"
    )
    args = sanitize_help(parser).parse_args()

    graphbase = args.graphbase

    # @RamRS: This might need some more work
    infiles = [graphbase, graphbase + '.tagset']
    if os.path.exists(graphbase + '.stoptags'):
        infiles.append(graphbase + '.stoptags')
    for _ in infiles:
        check_input_files(_, args.force)

    check_space(infiles, args.force)

    print('loading k-mer nodegraph %s' % graphbase, file=sys.stderr)
    graph = Nodegraph.load(graphbase)

    print('loading tagset %s.tagset...' % graphbase, file=sys.stderr)
    graph.load_tagset(graphbase + '.tagset')

    initial_stoptags = False    # @CTB regularize with make-initial
    if os.path.exists(graphbase + '.stoptags'):
        print('loading stoptags %s.stoptags' % graphbase, file=sys.stderr)
        graph.load_stop_tags(graphbase + '.stoptags')
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

    # create countgraph
    ksize = graph.ksize()
    counting = khmer_args.create_countgraph(args, ksize=ksize)

    # load & merge
    for index, subset_file in enumerate(pmap_files):
        print('<-', subset_file, file=sys.stderr)
        subset = SubsetPartition.load(subset_file, graph)

        print('** repartitioning subset... %s' % subset_file, file=sys.stderr)
        graph.repartition_largest_partition(counting,
                                            EXCURSION_DISTANCE,
                                            EXCURSION_KMER_THRESHOLD,
                                            EXCURSION_KMER_COUNT_THRESHOLD,
                                            subs=subset)

        print('** merging subset... %s' % subset_file, file=sys.stderr)
        graph.merge_subset(subset)

        print('** repartitioning, round 2... %s' %
              subset_file, file=sys.stderr)
        size = \
            graph.repartition_largest_partition(counting,
                                                EXCURSION_DISTANCE,
                                                EXCURSION_KMER_THRESHOLD,
                                                EXCURSION_KMER_COUNT_THRESHOLD)

        print('** repartitioned size:', size, file=sys.stderr)

        print('saving stoptags binary', file=sys.stderr)
        graph.save_stop_tags(graphbase + '.stoptags')
        os.rename(subset_file, subset_file + '.processed')
        print('(%d of %d)\n' % (index, len(pmap_files)), file=sys.stderr)

    print('done!', file=sys.stderr)


if __name__ == '__main__':
    main()
