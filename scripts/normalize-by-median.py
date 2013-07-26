#! /usr/bin/env python
"""
Eliminate reads with median k-mer abundance higher than
DESIRED_COVERAGE.  Output sequences will be placed in 'infile.keep'.

% python scripts/normalize-by-median.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.
"""

import sys
import screed
import os
import argparse

import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE
from khmer.threading_args import add_threading_args
from khmer import thread_utils
from khmer.thread_utils import ThreadedWriter, PairThreadedWriter

DEFAULT_DESIRED_COVERAGE = 10

def main():
    parser = build_construct_args()
    add_threading_args(parser)
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('-p', '--paired', action='store_true')
    parser.add_argument('-s', '--savehash', dest='savehash', default='')
    parser.add_argument('-l', '--loadhash', dest='loadhash',
                        default='')
    parser.add_argument('-R', '--report-to-file', dest='report_file',
                        type=argparse.FileType('w'))
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, \
                "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, \
            ' - min hashsize = %-5.2g \t(-x)' % args.min_hashsize
        print >>sys.stderr, ' - paired =	      %s \t\t(-p)' % args.paired
        print >>sys.stderr, ''
        print >>sys.stderr, \
            'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize)' \
            % (args.n_hashes * args.min_hashsize)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes
    DESIRED_COVERAGE = args.cutoff
    report_fp = args.report_file
    filenames = args.input_filenames
    n_threads = int(args.n_threads)
    if n_threads == 1:
        n_threads = 2

    if args.loadhash:
        print 'loading hashtable from', args.loadhash
        ht = khmer.load_counting_hash(args.loadhash)
    else:
        print 'making hashtable'
        ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    total = 0
    discarded = 0

    def single_filter_fn(name, seq):
        if len(seq) < K:
            return None

        seq = seq.replace('N', 'A')
        med, _, _ = ht.get_median_count(seq)

        if med < DESIRED_COVERAGE:
            ht.consume(seq)
            return seq

    def pair_filter_fn(n1, s1, n2, s2):
        if single_filter_fn(n1, s1) or \
           single_filter_fn(n2, s2):
            return s1, s2
        else:
            return None, None

    for input_filename in filenames:
        output_name = os.path.basename(input_filename) + '.keep'
        outfp = open(output_name, 'w')

        if args.paired:
            writer_class = PairThreadedWriter
            filter_fn = pair_filter_fn
        else:
            writer_class = ThreadedWriter
            filter_fn = single_filter_fn

        tw = writer_class(outfp).start()
        rparser = khmer.ReadParser(input_filename, n_threads - 1)
        threads = thread_utils.start_threads(n_threads - 1,
                                             target=tw.process_fn,
                                             args=(rparser, filter_fn))

        tw.join(threads)

        print 'DONE with', input_filename
        # end loop

    if args.savehash:
        print 'Saving hashfile through', input_filename
        print '...saving to', args.savehash
        ht.save(args.savehash)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the counting hash is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        print >>sys.stderr, "** Do not use these results!!"
        sys.exit(-1)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
