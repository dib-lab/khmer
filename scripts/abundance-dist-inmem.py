#! /usr/bin/env python
"""
Produce the k-mer abundance distribution for the given file, without
loading a prebuilt counting hash.

% python scripts/abundance-dist-inmem.py <data> <histout>

Use '-h' for parameter help.
"""
import sys
import khmer
import argparse
import os
import threading
from khmer.counting_args import build_construct_args, report_on_config
from khmer.threading_args import add_threading_args

def main():
    parser = build_construct_args(
        "Output k-mer abundance distribution (inmem version).")
    add_threading_args(parser)

    parser.add_argument('datafile')
    parser.add_argument('histout')

    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--savehash', dest='savehash', default='')

    args = parser.parse_args()
    report_on_config(args)

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    n_threads = int(args.n_threads)
    
    datafile = args.datafile
    histout = args.histout

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT, n_threads)
    ht.set_use_bigcount(args.bigcount)

    print 'building tracking ht'
    K = ht.ksize()
    sizes = ht.hashsizes()
    tracking = khmer._new_hashbits(K, sizes)

    print 'K:', K
    print 'HT sizes:', sizes
    print 'outputting to', histout

    # start loading
    rparser = khmer.ReadParser(datafile, n_threads)
    threads = []
    print 'consuming input', datafile
    for tnum in xrange(n_threads):
        t = \
            threading.Thread(
                target=ht.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(t)
        t.start()

    for t in threads:
        t.join()

    print 'preparing hist from %s...' % datafile
    z = ht.abundance_distribution(datafile, tracking)
    total = sum(z)

    if 0 == total:
        print >>sys.stderr, \
            "ERROR: abundance distribution is uniformly zero; " \
            "nothing to report."
        print >>sys.stderr, "\tPlease verify that the input files are valid."
        sys.exit(-1)

    fp = open(histout, 'w')

    sofar = 0
    for n, i in enumerate(z):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print >>fp, n, i, sofar, round(frac, 3)

        if sofar == total:
            break

    if args.savehash:
        print 'Saving hashfile', args.savehash
        print '...saving to', args.savehash
        ht.save(args.savehash)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
