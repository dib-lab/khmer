#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import division
from __future__ import print_function
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

    print('hashtable from', hashfile)
    ht = khmer.load_counting_hash(hashfile)

    hist = {}

    for i in range(65536):
        hist[i] = 0

    for n, record in enumerate(screed.open(seqfile)):
        if n > 0 and n % 100000 == 0:
            print('...', n)

        seq = record.sequence.replace('N', 'A')

        try:
            med, _, _ = ht.get_median_count(seq)
        except ValueError:
            continue

        hist[med] = hist[med] + 1

    histlist = list(hist.items())
    histlist.sort()

    maxk = max(hist.keys())
    sumk = sum(hist.values())

    sofar = 0
    for n, m in histlist:
        sofar += m
        percent = float(sofar) / sumk
        outfp.write('%d %d %d %.3f\n' % (n, m, sofar, percent))
    outfp.close()

if __name__ == '__main__':
    main()
