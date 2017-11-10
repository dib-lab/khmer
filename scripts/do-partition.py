#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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
# pylint: disable=missing-docstring,invalid-name
"""
Do all the partition steps in one script.

% do-partition.py <graphname> <reads1> [ <reads2> ... ]

Use '-h' for parameter help.
"""

import khmer
import sys
import threading
import os.path
import os
import textwrap
from khmer import khmer_args
from khmer.khmer_args import (build_nodegraph_args, report_on_config,
                              add_threading_args, sanitize_help)
import glob
from khmer.kfile import check_input_files, check_space
from oxli.partition import worker

# stdlib queue module was renamed on Python 3
try:
    import queue
except ImportError:
    import Queue as queue

DEFAULT_SUBSET_SIZE = int(1e5)
DEFAULT_N_THREADS = 4
DEFAULT_K = 32


def get_parser():
    epilog = """\
    Load in a set of sequences, partition them, merge the partitions, and
    annotate the original sequences files with the partition information.

    This script combines the functionality of
    :program:`load-graph.py`, :program:`partition-graph.py`,
    :program:`merge-partitions.py`, and :program:`annotate-partitions.py` into
    one script. This is convenient but should probably not be used for large
    data sets, because :program:`do-partition.py` doesn't provide save/resume
    functionality.

    Example::

        do-partition.py -k 20 example tests/test-data/random-20-a.fa
    """
    parser = build_nodegraph_args(
        descr='Load, partition, and annotate FAST[AQ] sequences',
        epilog=textwrap.dedent(epilog), citations=['graph'])
    add_threading_args(parser)
    parser.add_argument('-s', '--subset-size', default=DEFAULT_SUBSET_SIZE,
                        dest='subset_size', type=float,
                        help='Set subset size (usually 1e5-1e6 is good)')
    parser.add_argument('--no-big-traverse', dest='no_big_traverse',
                        action='store_true', default=False,
                        help='Truncate graph joins at big traversals')
    parser.add_argument('--keep-subsets', default=False, action='store_true',
                        help='Keep individual subsets')
    parser.add_argument('graphbase', help="base name for output files")
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+', help='input FAST[AQ] sequence filenames')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


# pylint: disable=too-many-branches
def main():  # pylint: disable=too-many-locals,too-many-statements
    args = sanitize_help(get_parser()).parse_args()

    report_on_config(args, graphtype='nodegraph')

    for infile in args.input_filenames:
        check_input_files(infile, args.force)

    check_space(args.input_filenames, args.force)

    print('Saving k-mer nodegraph to %s' %
          args.graphbase, file=sys.stderr)
    print('Loading kmers from sequences in %s' %
          repr(args.input_filenames), file=sys.stderr)
    print('--', file=sys.stderr)
    print('SUBSET SIZE', args.subset_size, file=sys.stderr)
    print('N THREADS', args.threads, file=sys.stderr)
    print('--', file=sys.stderr)

    # load-graph.py

    print('making nodegraph', file=sys.stderr)
    nodegraph = khmer_args.create_nodegraph(args)

    for _, filename in enumerate(args.input_filenames):
        print('consuming input', filename, file=sys.stderr)
        nodegraph.consume_seqfile_and_tag(filename)

    # 0.18 is ACTUAL MAX. Do not change.
    fp_rate = \
        khmer.calc_expected_collisions(
            nodegraph, args.force, max_false_pos=.15)
    print('fp rate estimated to be %1.3f' % fp_rate, file=sys.stderr)

    # partition-graph

    # do we want to exhaustively traverse the graph?
    stop_big_traversals = args.no_big_traverse
    if stop_big_traversals:
        print('** This script brakes for lumps: ',
              'stop_big_traversals is true.', file=sys.stderr)
    else:
        print('** Traverse all the things:',
              ' stop_big_traversals is false.', file=sys.stderr)

    #
    # now, partition!
    #

    # divide the tags up into subsets
    divvy = nodegraph.divide_tags_into_subsets(int(args.subset_size))
    divvy = list(divvy)
    n_subsets = len(divvy)
    divvy.append(0)

    # build a queue of tasks:
    worker_q = queue.Queue()

    # break up the subsets into a list of worker tasks
    for _ in range(0, n_subsets):
        start = divvy[_]
        end = divvy[_ + 1]
        worker_q.put((nodegraph, _, start, end))

    print('enqueued %d subset tasks' % n_subsets, file=sys.stderr)
    open('%s.info' % args.graphbase, 'w').write('%d subsets total\n'
                                                % (n_subsets))

    if n_subsets < args.threads:
        args.threads = n_subsets

    # start threads!
    print('starting %d threads' % args.threads, file=sys.stderr)
    print('---', file=sys.stderr)

    threads = []
    for _ in range(args.threads):
        cur_thread = threading.Thread(target=worker,
                                      args=(worker_q, args.graphbase,
                                            stop_big_traversals))
        threads.append(cur_thread)
        cur_thread.start()

    print('done starting threads', file=sys.stderr)

    # wait for threads
    for _ in threads:
        _.join()

    print('---', file=sys.stderr)
    print('done making subsets! see %s.subset.*.pmap' %
          (args.graphbase,), file=sys.stderr)

    # merge-partitions

    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print('loading %d pmap files (first one: %s)' %
          (len(pmap_files), pmap_files[0]), file=sys.stderr)

    nodegraph = khmer.Nodegraph(args.ksize, 1, 1)

    for pmap_file in pmap_files:
        print('merging', pmap_file, file=sys.stderr)
        nodegraph.merge_subset_from_disk(pmap_file)

    if not args.keep_subsets:
        print('removing pmap files', file=sys.stderr)
        for pmap_file in pmap_files:
            os.unlink(pmap_file)

    # annotate-partitions

    for infile in args.input_filenames:
        print('outputting partitions for', infile, file=sys.stderr)
        outfile = os.path.basename(infile) + '.part'
        part_count = nodegraph.output_partitions(infile, outfile)
        print('output %d partitions for %s' % (
            part_count, infile), file=sys.stderr)
        print('partitions are in', outfile, file=sys.stderr)


if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
