#! /usr/bin/env python
"""
Build a graph from the given sequences, save in <htname>.

% python scripts/load-graph.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import threading
import Queue
import os.path

import khmer
from khmer.hashbits_args import build_construct_args
from khmer.counting_args import report_on_config
from khmer.threading_args import add_threading_args


def main():
    parser = build_construct_args()
    add_threading_args(parser)
    parser.add_argument('--no-build-tagset', '-n', default=False,
                        action='store_true', dest='no_build_tagset',
                        help='Do NOT construct tagset while loading sequences')
    parser.add_argument('output_filename')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()
    report_on_config(args)

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    base = args.output_filename
    filenames = args.input_filenames
    n_threads = int(args.n_threads)

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)
    if args.no_build_tagset:
        print 'We WILL NOT build the tagset.'
    else:
        print 'We WILL build the tagset (for partitioning/traversal).'

    ###

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    if args.no_build_tagset:
        target_method = ht.consume_fasta_with_reads_parser
    else:
        target_method = ht.consume_fasta_and_tag_with_reads_parser

    for n, filename in enumerate(filenames):
        
        rparser = khmer.ReadParser(filename, n_threads)
        threads = []
        print 'consuming input', filename
        for tnum in xrange(n_threads):
            t = threading.Thread(target=target_method, args=(rparser, ))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

    print 'saving hashtable in', base + '.ht'
    ht.save(base + '.ht')

    if not args.no_build_tagset:
        print 'saving tagset in', base + '.tagset'
        ht.save_tagset(base + '.tagset')

    info_fp = open(base + '.info', 'w')
    info_fp.write('%d unique k-mers' % ht.n_unique_kmers())

    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be %1.3f' % fp_rate
    if fp_rate > 0.15:          # 0.18 is ACTUAL MAX. Do not change.
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the graph structure is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        sys.exit(-1)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
