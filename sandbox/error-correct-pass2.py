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

% python scripts/error-correct-pass2 <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import screed
import os
import khmer
import argparse


###

DEFAULT_CUTOFF = 2

def output_single(read, new_sequence):
    name = read.name
    sequence = new_sequence

    quality = None
    if hasattr(read, 'quality'):
        quality = read.quality[:len(sequence)]
        sequence = sequence[:len(quality)] # in cases where sequence _lengthened_

    if quality:
        assert len(sequence) == len(quality), (sequence, quality)
        return "@%s\n%s\n+\n%s\n" % (name, sequence, quality)
    else:
        return ">%s\n%s\n" % (name, sequence)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--trusted-cov", dest="trusted_cov", type=int,
                        default=DEFAULT_CUTOFF)
    parser.add_argument("--theta", dest="bits_theta", type=float, default=1.0)
    parser.add_argument('-o', '--output', dest='output_file',
                        help="output file for histogram; defaults to "
                             "<first filename>.errhist in cwd.",
                        type=argparse.FileType('w'), default=None)

    parser.add_argument('counts_table')
    parser.add_argument('readfile')
    
    args = parser.parse_args()

    print 'loading counts'
    ht = khmer.load_counting_hash(args.counts_table)

    aligner = khmer.ReadAligner(ht,
                                args.trusted_cov,
                                args.bits_theta)

    print "trusted:", args.trusted_cov

    corrfp = args.output_file
    if not corrfp:
        outfile = os.path.basename(args.readfile) + '.corr'
        corrfp = open(outfile, 'w')

    n_corrected = 0
    for n, read in enumerate(screed.open(args.readfile)):
        if n % 10000 == 0:
            print >>sys.stderr, '...', n, n_corrected
        seq = read.sequence.replace('N', 'A')

        # build the alignment...
        score, graph_alignment, read_alignment, truncated = \
               aligner.align(seq)
        
        if not truncated:
            graph_seq = graph_alignment.replace("-", "")
            if graph_seq != seq:
                n_corrected += 1

            seq = graph_seq

        corrfp.write(output_single(read, seq))

if __name__ == '__main__':
    main()
