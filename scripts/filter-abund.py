#! /usr/bin/env python
"""
Trim sequences at k-mers of the given abundance, based on the given counting
hash table.  Output sequences will be placed in 'infile.abundfilt'.

% python scripts/filter-abund.py <counting.kh> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer import threading_args as targs
from khmer.counting_args import build_counting_multifile_args

###

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2

def main():
    parser = build_counting_multifile_args()
    targs.add_threading_args(parser)
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
    n_threads = int(args.n_threads)

    print 'file with ht: %s' % counting_ht

    print 'loading hashtable'
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    print "K:", K

    ### the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        if args.variable_coverage: # only trim when sequence has high enough C
            med, _, _ = ht.get_median_count(seq)
            if med < args.normalize_to:
                return name, seq

        trim_seq, trim_at = ht.trim_on_abundance(seq, args.cutoff)

        if trim_at >= K:
            return name, trim_seq

        return None, None

    ### the filtering loop
    for infile in infiles:
        print 'filtering', infile
        outfile = os.path.basename(infile) + '.abundfilt'
        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn, n_workers=n_threads)
        tsp.start(verbose_loader(infile), outfp)

        print 'output in', outfile

if __name__ == '__main__':
    main()
