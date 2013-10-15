#! /usr/bin/env python
import khmer
import sys
import argparse
import os, os.path
import screed

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('hashname')
    parser.add_argument('datafiles', nargs='+')

    args = parser.parse_args()
    hashfile = args.hashname
    datafiles = args.datafiles

    print 'loading counting hash'
    ht = khmer.load_counting_hash(hashfile)
    print 'loaded.'

    for datafile in datafiles:
        print 'annotating', datafile

        outfile = os.path.basename(datafile) + '.kannot'
        outfp = open(outfile, 'w')

        for n, record in enumerate(screed.open(datafile)):
            if n % 1000 == 0:
                print '...', n
            med, _, _ = ht.get_median_count(record.sequence)
            outfp.write('>%s kmed=%d\n%s\n' % (record.name, med,
                                               record.sequence))

if __name__ == '__main__':
    main()
