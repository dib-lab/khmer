#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import khmer
import argparse
import os
import screed


def main():
    parser = argparse.ArgumentParser(
        description="Output k-mer abundance distribution.")

    parser.add_argument('hashname')
    parser.add_argument('seqfile')
    parser.add_argument('histout')

    args = parser.parse_args()
    hashfile = args.hashname
    seqfile = args.seqfile
    histout = args.histout

    outfp = open(histout, 'w')

    print 'hashtable from', hashfile
    ht = khmer.load_counting_hash(hashfile)

    hist = {}

    for n, record in enumerate(screed.open(seqfile)):
        if n > 0 and n % 100000 == 0:
            print '...', n

        seq = record.sequence.replace('N', 'A')
        med, _, _ = ht.get_median_count(seq)

        hist[med] = hist.get(med, 0) + 1

    maxk = max(hist.keys())

    for i in range(maxk + 1):
        outfp.write('%d %d\n' % (i, hist.get(i, 0)))
    outfp.close()

if __name__ == '__main__':
    main()
