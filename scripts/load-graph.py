#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Build a graph from the given sequences, save in <htname>.

% python scripts/load-graph.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import threading

import khmer
from khmer.khmer_args import build_hashbits_args
from khmer.khmer_args import (report_on_config, info)
from khmer.threading_args import add_threading_args
from khmer.file import check_file_status, check_space
from khmer.file import check_space_for_hashtable


def get_parser():
    parser = build_hashbits_args(descr="Load sequences into the compressible "
                                 "graph format plus optional tagset.")
    add_threading_args(parser)
    parser.add_argument('--no-build-tagset', '-n', default=False,
                        action='store_true', dest='no_build_tagset',
                        help='Do NOT construct tagset while loading sequences')
    parser.add_argument('output_filename')
    parser.add_argument('input_filenames', nargs='+')
    return parser


def main():
    info('load-graph.py', ['graph'])
    args = get_parser().parse_args()
    report_on_config(args, hashtype='hashbits')

    base = args.output_filename
    filenames = args.input_filenames
    n_threads = int(args.n_threads)

    for _ in args.input_filenames:
        check_file_status(_)

    check_space(args.input_filenames)
    check_space_for_hashtable(args.ksize * args.min_tablesize)

    print 'Saving k-mer presence table to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)
    if args.no_build_tagset:
        print 'We WILL NOT build the tagset.'
    else:
        print 'We WILL build the tagset (for partitioning/traversal).'

    print 'making k-mer presence table'
    htable = khmer.new_hashbits(args.ksize, args.min_tablesize, args.n_tables)

    if args.no_build_tagset:
        target_method = htable.consume_fasta_with_reads_parser
    else:
        target_method = htable.consume_fasta_and_tag_with_reads_parser

    config = khmer.get_config()
    config.set_reads_input_buffer_size(n_threads * 64 * 1024)

    for _, filename in enumerate(filenames):

        rparser = khmer.ReadParser(filename, n_threads)
        threads = []
        print 'consuming input', filename
        for _ in xrange(n_threads):
            cur_thrd = threading.Thread(target=target_method, args=(rparser, ))
            threads.append(cur_thrd)
            cur_thrd.start()

        for thread in threads:
            thread.join()

    print 'saving k-mer presence table in', base + '.pt'
    htable.save(base + '.pt')

    if not args.no_build_tagset:
        print 'saving tagset in', base + '.tagset'
        htable.save_tagset(base + '.tagset')

    info_fp = open(base + '.info', 'w')
    info_fp.write('%d unique k-mers' % htable.n_unique_kmers())

    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be %1.3f' % fp_rate
    if fp_rate > 0.15:          # 0.18 is ACTUAL MAX. Do not change.
        print >> sys.stderr, "**"
        print >> sys.stderr, ("** ERROR: the graph structure is too small for "
                              "this data set.  Increase table size/# tables.")
        print >> sys.stderr, "**"
        sys.exit(1)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
