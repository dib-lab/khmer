#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
"""


import argparse
import math
import os
import sys
import textwrap

import khmer
from khmer.khmer_args import DEFAULT_K, info, ComboFormatter
from khmer import __version__
from screed.fasta import fasta_iter


def get_parser():
    descr = ""
    epilog = ("")
    parser = argparse.ArgumentParser(
        description=descr, epilog=textwrap.dedent(epilog),
        formatter_class=ComboFormatter)

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)

    parser.add_argument('--version', action='version',
                        version='khmer {v}'.format(v=__version__))
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')

    parser.add_argument('--ksize', '-k', type=int, default=env_ksize,
                        help='k-mer size to use')

    parser.add_argument('--error-rate', '-e', type=float, default=0.01,
                        help='Acceptable error rate')

    parser.add_argument('-R', '--report',
                        metavar='filename', type=argparse.FileType('w'))

    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')


    return parser


def main():
    info('kmer-intersection.py')
    args = get_parser().parse_args()
    report_fp = args.report
    input_filename = None

    total_hll = khmer.HLLCounter(args.error_rate, args.ksize)
    curve = []

    first = khmer.HLLCounter(args.error_rate, args.ksize)
    first.consume_fasta(args.input_filenames[0])
    total_hll.merge(first)

    second = khmer.HLLCounter(args.error_rate, args.ksize)
    for n, record in enumerate(fasta_iter(open(args.input_filenames[1]))):
        second.consume_string(record['sequence'])
        interval = int(math.log(n + 1, 1.1))
        if n < 100 or n % interval == 0:
            total_hll.merge(second)
            curve.append((n, len(first) + len(second) - len(total_hll)))

    overlap = len(first) + len(second) - len(total_hll)
    print '# of unique k-mers in dataset 1:', len(first)
    print '# of unique k-mers in dataset 2:', len(second)
    print '# of overlap unique k-mers:', overlap

    total_reads = curve[-1][0]
    interval = total_reads / 100
    with open('curve', 'w') as f:
        for c in curve:
            if c[0] % interval == 0:  # TODO: this is WRONG!
                f.write("%d %d\n" % (c[0], c[1]))


if __name__ == "__main__":
    main()
