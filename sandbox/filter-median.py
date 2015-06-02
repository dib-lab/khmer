#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Accept or discard sequences XXX, based on the given counting
hash table.  Output sequences will be placed in 'infile.medfilt'.

% python sandbox/filter-median.py <counting.ct> <data1> [ <data2> <...> ]

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
    parser = build_counting_args()
    parser.add_argument('--coverage', '-C', dest='coverage',
                        default=DEFAULT_COVERAGE, type=int)
    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print('file with ht: %s' % counting_ht)

    print('loading hashtable')
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    print("K:", K)

    # the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']

        med, avg, dev = ht.get_median_count(seq)

        if random.randint(1, med) > args.coverage:
            return None, None

        return name, seq

    # the filtering loop
    for infile in infiles:
        print('filtering', infile)
        outfile = os.path.basename(infile) + '.medfilt'
        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), outfp)

        print('output in', outfile)

if __name__ == '__main__':
    main()
