#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring,no-member
"""
Estimate number of unique k-mers, with precision <= ERROR_RATE.

% python scripts/unique-kmers.py [ -k <k size> ] [ -e <ERROR_RATE> ] <data1>
<data2> ...

Use '-h' for parameter help.
"""
from __future__ import print_function


import argparse
import os
import sys
import textwrap

import khmer
from khmer.khmer_args import (DEFAULT_K, info, ComboFormatter,
                              _VersionStdErrAction)
from khmer.utils import write_record
from khmer.khmer_args import graphsize_args_report
from khmer import __version__
import screed


def get_parser():
    descr = "Estimate number of unique k-mers, with precision <= ERROR_RATE."
    epilog = ("""
    A HyperLogLog counter is used to do cardinality estimation. Since this counter
    is based on a tradeoff between precision and memory consumption,
    the :option:`-e`/:option:`--error-rate` can be used to control how much
    memory will be used. In practice the memory footprint is small even
    at low error rates (< 0.01).

    :option:`-k`/:option:`--ksize` should be set to the desired k-mer size.

    Informational output is sent to STDERR, but a report file can be generated
    with :option:`-R`/:option:`--report`.

    :option:`--stream-out` will write the sequences taken in to STDOUT.
    This is useful for workflows: count unique kmers in a stream, then do
    digital normalization.

    :option:`--diagnostics` will provide detailed options for tablesize
    and memory limitations for various false positive rates. This is useful for
    configuring other khmer scripts. This will be written to STDERR.

    Example::

        unique-kmers.py -k 17 tests/test-data/test-abund-read{,-2,-3}.fa

    Example::

        unique-kmers.py -k 17 --diagnostics tests/test-data/test-abund-read.fa

    Example::

        unique-kmers.py --stream-out -k 17 tests/test-data/test-reads.fa | \\
        normalize-by-median.py -k 17 -o normalized /dev/stdin

    Example::

        unique-kmers.py -R unique_count -k 30 \\
        tests/test-data/test-abund-read-paired.fa""")  # noqa
    parser = argparse.ArgumentParser(
        description=descr, epilog=textwrap.dedent(epilog),
        formatter_class=ComboFormatter)

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)

    parser.add_argument('--version', action=_VersionStdErrAction,
                        version='khmer {v}'.format(v=__version__))

    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')

    parser.add_argument('--ksize', '-k', type=int, default=env_ksize,
                        help='k-mer size to use')

    parser.add_argument('--error-rate', '-e', type=float, default=0.01,
                        help='Acceptable error rate')

    parser.add_argument('--report', '-R',
                        metavar='filename', type=argparse.FileType('w'),
                        help='generate informational report and write to'
                        ' filename')

    parser.add_argument('--stream-out', '-S', default=False,
                        action='store_true',
                        help='write input sequences to STDOUT')

    parser.add_argument('--diagnostics', default=False, action='store_true',
                        help='print out recommended tablesize arguments and '
                             'restrictions')

    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename(s).', nargs='+')

    return parser


def main():
    info('unique-kmers.py', ['SeqAn', 'hll'])
    args = get_parser().parse_args()

    total_hll = khmer.HLLCounter(args.error_rate, args.ksize)

    report_fp = args.report
    input_filename = None
    for index, input_filename in enumerate(args.input_filenames):
        hllcpp = khmer.HLLCounter(args.error_rate, args.ksize)
        for record in screed.open(input_filename):
            seq = record.sequence.upper().replace('N', 'A')
            hllcpp.consume_string(seq)
            if args.stream_out:
                write_record(record, sys.stdout)

        cardinality = hllcpp.estimate_cardinality()
        print('Estimated number of unique {0}-mers in {1}: {2}'.format(
              args.ksize, input_filename, cardinality),
              file=sys.stderr)

        if report_fp:
            print(cardinality, args.ksize, '(total)', file=report_fp)
            report_fp.flush()
        total_hll.merge(hllcpp)

    cardinality = total_hll.estimate_cardinality()
    print('Total estimated number of unique {0}-mers: {1}'.format(
          args.ksize, cardinality),
          file=sys.stderr)

    to_print = graphsize_args_report(cardinality, args.error_rate)
    if args.diagnostics:
        print(to_print, file=sys.stderr)

    if report_fp:
        print(cardinality, args.ksize, 'total', file=report_fp)
        print(to_print, file=report_fp)
        report_fp.flush()

if __name__ == "__main__":
    main()
