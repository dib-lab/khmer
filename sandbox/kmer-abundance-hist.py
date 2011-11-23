#! /usr/bin/env python
"""
Produce the k-mer abundance distribution for the given file.

% python scripts/kmer-abundance-hist.py [ -z ] <htname> <data> <histout>

Use '-h' for parameter help.
"""
import sys, khmer
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Output k-mer abundance distribution.")
    
    parser.add_argument('hashname')
    parser.add_argument('datafile')
    parser.add_argument('histout')

    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')

    args = parser.parse_args()
    hashfile = args.hashname
    datafile = args.datafile
    histout = args.histout

    print 'hashtable from', hashfile
    ht = khmer.load_counting_hash(hashfile)

    K = ht.ksize()
    sizes = ht.hashsizes()
    tracking = khmer._new_hashbits(K, sizes)

    print 'K:', K
    print 'HT sizes:', sizes
    print 'outputting to', histout

    if os.path.exists(histout):
        if not args.squash_output:
            print >>sys.stderr, 'ERROR: %s exists; not squashing.' % histout
            sys.exit(-1)
        
        print '** squashing existing file %s' % histout

    print 'preparing hist...'
    z = ht.abundance_distribution(datafile, tracking)
    total = sum(z[1:])
        
    fp = open(histout, 'w')
    sofar = 0
    for n, i in enumerate(z[1:]):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print >>fp, n + 1, i, sofar, round(frac, 3)

if __name__ == '__main__':
    main()
