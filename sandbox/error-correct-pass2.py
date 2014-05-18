#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Error correct reads based on a counting hash from a diginorm step.
Output sequences will be put in @@@.

% python scripts/error-correct-pass2 <counting.kh> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

from khmer.counting_args import build_counting_multifile_args

###

DEFAULT_COVERAGE = 20
DEFAULT_MAX_ERROR_REGION = 40


def main():
    parser = build_counting_multifile_args()
    parser.add_argument('--cutoff', '-C', dest='coverage',
                        default=DEFAULT_COVERAGE, type=int,
                        help="Diginorm coverage.")
    parser.add_argument('--max-error-region', '-M', dest='max_error_region',
                        default=DEFAULT_MAX_ERROR_REGION, type=int,
                        help="Max length of error region allowed")
    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print 'file with ht: %s' % counting_ht

    print 'loading hashtable'
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()
    C = args.coverage
    max_error_region = args.max_error_region

    print "K:", K
    print "C:", C
    print "max error region:", max_error_region

    # the filtering function.
    def process_fn(record):
        # read_aligner is probably not threadsafe?
        aligner = khmer.new_readaligner(ht, 1, C, max_error_region)

        name = record['name']
        seq = record['sequence']

        seq = seq.replace('N', 'A')

        grXreAlign, reXgrAlign = aligner.align(seq)

        if len(reXgrAlign) > 0:
            graph_seq = grXreAlign.replace('-', '')
            seq = graph_seq

        return name, seq

    # the filtering loop
    for infile in infiles:
        print 'filtering', infile
        outfile = os.path.basename(infile) + '.corr'
        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), outfp)

        print 'output in', outfile

if __name__ == '__main__':
    main()
