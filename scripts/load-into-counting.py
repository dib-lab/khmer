#! /usr/bin/env python
"""
Build a counting Bloom filter from the given sequences, save in <htname>.

% python scripts/load-into-counting.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.

@CTB enable/disable bigcount via command-line parameters.
"""

import sys
import screed
import threading
import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE
from khmer.threading_args import add_threading_args

###


def main():
    parser = build_construct_args()
    add_threading_args(parser)
    parser.add_argument('output_filename')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print>>sys.stderr, "** WARNING: hashsize is default!  You absodefly want to increase this!\n** Please read the docs!"

        print>>sys.stderr, '\nPARAMETERS:'
        print>>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print>>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print>>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % args.min_hashsize
        print>>sys.stderr, ''
        print>>sys.stderr, 'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize)' % (
            args.n_hashes * args.min_hashsize)
        print>>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    base = args.output_filename
    filenames = args.input_filenames
    n_threads = int(args.n_threads)

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)

    ###

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT, n_threads)
    ht.set_use_bigcount(True)

    for n, filename in enumerate(filenames):

        rparser = khmer.ReadParser(filename, n_threads)
        threads = []
        print 'consuming input', filename
        for tnum in xrange(n_threads):
            t = \
                threading.Thread(
                    target=ht.consume_fasta_with_reads_parser,
                    args=(rparser, )
                )
            threads.append(t)
            t.start()
            # ht.consume_fasta(filename)

        for t in threads:
            t.join()

        if n > 0 and n % 10 == 0:
            print 'mid-save', base
            ht.save(base)
            open(base + '.info', 'w').write('through %s' % filename)

    print 'saving', base
    ht.save(base)

    info_fp = open(base + '.info', 'w')
    info_fp.write('through end: %s\n' % filename)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print 'fp rate estimated to be %1.3f' % fp_rate
    print >>info_fp, 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the counting hash is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        sys.exit(-1)

    print 'DONE.'

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
