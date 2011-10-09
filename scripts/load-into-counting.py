# /usr/bin/env python
"""
Build a counting Bloom filter from the given sequences, save in <htname>.

% python scripts/load-into-counting.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys, screed
import khmer
from khmer.counting_args import build_construct_args, DEFAULT_MIN_HASHSIZE

###

def main():
    parser = build_construct_args()
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
        print>>sys.stderr, 'Estimated memory usage is %.2g bytes (n_hashes x min_hashsize)' % (args.n_hashes * args.min_hashsize)
        print>>sys.stderr, '-'*8


    K=args.ksize
    HT_SIZE=args.min_hashsize
    N_HT=args.n_hashes

    base = args.output_filename
    filenames = args.input_filenames

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)

    ###
    
    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)
    ht.set_use_bigcount(True)

    for n, filename in enumerate(filenames):
       print 'consuming input', filename
       ht.consume_fasta(filename)

       if n > 0 and n % 10 == 0:
           print 'mid-save', base
           ht.save(base)
           open(base + '.info', 'w').write('through %s' % filename)

    print 'saving', base
    ht.save(base)
    open(base + '.info', 'w').write('through end: %s' % filename)

if __name__ == '__main__':
    main()
