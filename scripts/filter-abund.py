#! /usr/bin/env python
"""
Trim sequences at k-mers of the given abundance, based on the given counting
hash table.  Output sequences will be placed in 'infile.abundfilt'.

% python scripts/filter-abund.py <counting.kh> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import os
import threading

import khmer

from khmer.counting_args import build_counting_multifile_args
from khmer.threading_args import add_threading_args
from khmer import thread_utils
from khmer.thread_utils import ThreadedProcessor

###

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2

def main():
    parser = build_counting_multifile_args()
    add_threading_args(parser)
    parser.add_argument('--cutoff', '-C', dest='cutoff',
                        default=DEFAULT_CUTOFF, type=int,
                        help="Trim at k-mers below this abundance.")

    parser.add_argument('-V', '--variable-coverage', action='store_true',
                        dest='variable_coverage', default=False)
    parser.add_argument('--normalize-to', '-Z', type=int, dest='normalize_to',
                        help='base variable-coverage cutoff on this median k-mer abundance',
                        default=DEFAULT_NORMALIZE_LIMIT)

    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames
    n_threads = max(int(args.n_threads), 2) # one reader, one writer

    print 'file with ht: %s' % counting_ht
    print 'n_threads:', n_threads

    print 'loading hashtable'
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    print "K:", K
    config = khmer.get_config()
    config.set_reads_input_buffer_size(n_threads * 64 * 1024)

    ### the filtering function - return sequence, rest is taken care of.
    def filter_fn(name, seq):
        if 'N' in seq or len(seq) < K:
            return None

        # if var cov, only trim when sequence has high enough coverage
        if args.variable_coverage:
            med, _, _ = ht.get_median_count(seq)
            if med < args.normalize_to:
                return seq

        # high-coverage sequence or no var cov flag -- trim!
        trim_seq, trim_at = ht.trim_on_abundance(seq, args.cutoff)

        # is trimmed sequence long enough to keep?
        if trim_at >= K:
            return trim_seq

        # else, eliminate sequence
        return None
        
    ### the filtering loop
    for n, filename in enumerate(infiles):
        print 'filtering', filename
        outfile = os.path.basename(filename) + '.abundfilt'
        
        print 'output in', outfile
        outfp = open(outfile, 'w')

        # create and start a threaded writer
        tw = ThreadedProcessor(outfp).start()

        # create multithreaded readparser
        rparser = khmer.ReadParser(filename, n_threads - 1)
        threads = thread_utils.start_threads(n_threads - 1,
                                             target=tw.process_fn,
                                             args=(rparser, filter_fn))

        # wait for threads to finish & flush out any remaining records.
        if not tw.join(threads):
            print >>sys.stderr, "** failure.  See message above."
            sys.exit(-1)

if __name__ == '__main__':
    main()
