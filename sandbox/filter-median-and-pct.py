#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Accept or discard sequences XXX, based on the given counting
hash table.  Output sequences will be placed in 'infile.medpctfilt'.

% python sandbox/filter-median-and-pct.py <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
from __future__ import print_function
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

from khmer.khmer_args import build_counting_args

import random

###

DEFAULT_COVERAGE = 20


def main():
    parser = build_counting_multifile_args()
    parser.add_argument('--coverage', '-C', dest='coverage',
                        default=DEFAULT_COVERAGE, type=int)
    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print('file with ht: %s' % counting_ht)

    print('loading hashtable')
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    xxxfp = None

    print("K:", K)

    # the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']

        med, avg, dev = ht.get_median_count(seq)
        pct = dev / avg * 100

        xxxfp.write('%s %s %s %s %s\n' % (med, avg, dev, pct, name))

        if random.randint(1, med) > args.coverage or pct > 100:
            return None, None

        return name, seq

    # the filtering loop
    for infile in infiles:
        print('filtering', infile)
        xxxfp = open(os.path.basename(infile) + '.medpctfilt.stats', 'w')
        outfile = os.path.basename(infile) + '.medpctfilt'
        outfp = open(outfile, 'w')

        for n, record in enumerate(screed.open(infile)):
            if n % 100000 == 0:
                print('...', n)

            name, seq = process_fn(record)
            if name and seq:
                print('>%s\n%s' % (name, seq), file=outfp)

        print('output in', outfile)

if __name__ == '__main__':
    main()
