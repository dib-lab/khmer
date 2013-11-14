#! /usr/bin/env python
"""
Error correct reads based on a counting hash from a diginorm step.
Output sequences will be put in @@@.

% python scripts/error-correct-pass2 <counting.kh> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import screed
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

from khmer.counting_args import build_counting_multifile_args

###

DEFAULT_COVERAGE = 20
DEFAULT_MAX_ERROR_REGION = 40


def main():
    parser = build_counting_multifile_args()
    parser.add_argument("--trusted-cov", dest="trusted_cov", type=int, default=2)
    parser.add_argument("--theta", type=float, default=1.0)
    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print 'file with ht: %s' % counting_ht

    print 'loading hashtable'
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    aligner = khmer.new_readaligner(ht, args.trusted_cov, args.theta) # counting hash, trusted kmer coverage cutoff, bits theta (threshold value for terminating unproductive alignemnts)
    
    ### the filtering loop
    for infile in infiles:
        print 'aligning', infile
        for n, record in enumerate(screed.open(infile)):

            name = record['name']
            seq = record['sequence'].upper()
            print name
            print seq

            score, graph_alignment, read_alignment, truncated = aligner.align(seq)
            print score
            print graph_alignment
            print read_alignment
            print truncated

if __name__ == '__main__':
    main()
