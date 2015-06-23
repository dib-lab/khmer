#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Estimate optimal arguments using HLL couter.

% python sandbox/optimal_args_HLL.py [ -k <k size> ] [ -e <ERROR_RATE> ] <data1> <data2> ...
Use '-h' for parameter help.
"""

from __future__ import print_function

import argparse
import os
import sys
import textwrap
import math
import khmer
from khmer.khmer_args import DEFAULT_K, info, ComboFormatter
from khmer import __version__
from oxli.functions import optimal_args_output_gen as output_gen

def get_parser():
    descr = "Estimate optimal arguments using HLL couter., with precision <= ERROR_RATE."
    epilog = ("""
    A HyperLogLog counter is used to do cardinality estimation. Since this counter
    is based on a tradeoff between precision and memory consumption,
    :option:`-e`/:option:`--error-rate` can be used to control how much
    memory will be used. In practice the memory footprint is small even
    at low error rates (< 0.01).
    :option:`-k`/:option:`--ksize` should be set to the desired k-mer size.
    Output is sent to STDOUT, but a report file can be generated with
    :option:`-R`/:option:`--report`.
    Example::
        optimal_args_HLL.py -k 17 tests/test-data/test-abund-read{,-2,-3}.fa
    Example::
""" "    optimal_args_HLL.py.py -R optimal_args -k 30 tests/test-data/test-abund-read-paired.fa")  # noqa
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
    info('optimal_args_HLL.py', ['SeqAn'])
    args = get_parser().parse_args()

    hllcpp = khmer.HLLCounter(args.error_rate, args.ksize)

    report_fp = args.report
    fp_rate = args.error_rate
    input_filename = None
    for index, input_filename in enumerate(args.input_filenames):
        hllcpp.consume_fasta(input_filename)

    cardinality = hllcpp.estimate_cardinality()
    print('Estimated number of unique k-mers: {0}'.format(cardinality),
          file=sys.stderr)
        
    to_print = output_gen(cardinality,fp_rate)
    
    print(to_print)
    if report_fp:
        print(to_print, file=report_fp)
        report_fp.flush()


if __name__ == "__main__":
    main()
