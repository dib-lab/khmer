#! /usr/bin/env python
"""
Use a set of query reads to sweep out overlapping reads from multiple files.

% python scripts/sweep-reads3.py <query1> [ <query2> ... ] <search reads>

Results end up in <search reads>.sweep3.

Use '-h' for parameter help.
"""

import sys, khmer
import os.path
import screed
from khmer.hashbits_args import build_construct_args, DEFAULT_MIN_HASHSIZE

def main():
    parser = build_construct_args()
    parser.add_argument('input_filenames', nargs='+')
    parser.add_argument('read_filename')

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

    inputlist = args.input_filenames
    readsfile = args.read_filename

    query_list = []
    for n, inp_name in enumerate(inputlist):
        # create a hashbits data structure
        ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

        outfile = os.path.basename(inp_name) + '.sweep3'
        outfp = open(outfile, 'w')
        query_list.append((ht, outfp))

    for n, inp_name in enumerate(inputlist):
        ht = query_list[n][0]

        # load contigs, connect into N partitions
        print 'loading input reads from', inp_name
        ht.consume_fasta(inp_name)

    print 'starting sweep.'

    n = 0
    m = 0
    for n, record in enumerate(screed.open(readsfile)):
        if n % 10000 == 0:
            print '...', n, m

        for ht, outfp in query_list:
            count = ht.get_median_count(record.sequence)[0]
            if count:
                outfp.write('>%s\n%s\n' % (record.name, record.sequence))

if __name__ == '__main__':
    main()
