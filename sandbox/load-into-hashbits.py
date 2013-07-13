#! /usr/bin/env python
"""
Build a graph Bloom filter from the given sequences, save in <htname>.

% python sandbox/load-into-hashbits.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import screed
import os
import khmer
from khmer.counting_args import build_construct_args, report_on_config

def main():
    parser = build_construct_args()
    add_threading_args(parser)
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

    ###

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

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

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
