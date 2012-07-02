#! /usr/bin/env python
"""
Build a graph from the given sequences, save in <htname>.

% python scripts/load-graph.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import khmer, sys
import threading
import Queue
import gc
import os.path

import sys, screed
import khmer
from khmer.hashbits_args import build_construct_args, DEFAULT_MIN_HASHSIZE

def main():
    parser = build_construct_args()
    parser.add_argument('--no-build-tagset', '-n', default=False,
                        action='store_true', dest='no_build_tagset',
                        help='Do NOT construct tagset while loading sequences')
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
        print>>sys.stderr, 'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize / 8)' % (args.n_hashes * args.min_hashsize / 8.)
        print>>sys.stderr, '-'*8

    K=args.ksize
    HT_SIZE=args.min_hashsize
    N_HT=args.n_hashes

    base = args.output_filename
    filenames = args.input_filenames

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)
    if args.no_build_tagset:
        print 'We WILL NOT build the tagset.'
    else:
        print 'We WILL build the tagset (for partitioning/traversal).'

    ###
    
    print 'making hashtable'
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    for n, filename in enumerate(filenames):
       print 'consuming input', filename
       if args.no_build_tagset:
           ht.consume_fasta(filename)
       else:
           ht.consume_fasta_and_tag(filename)

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
