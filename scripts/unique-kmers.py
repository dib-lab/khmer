#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=invalid-name,missing-docstring,no-member
"""
Estimate number of unique k-mers, with precision <= ERROR_RATE.

% python scripts/unique-kmers.py [ -k <k size> ] [ -e <ERROR_RATE> ] <data1>
<data2> ...

Use '-h' for parameter help.
"""

import argparse
import os
import sys
import textwrap

import khmer
from khmer.khmer_args import DEFAULT_K, sanitize_help, KhmerArgumentParser
from khmer.khmer_args import graphsize_args_report


def get_parser():
    descr = "Estimate number of unique k-mers, with precision <= ERROR_RATE."
    epilog = """\
    A HyperLogLog counter is used to do cardinality estimation. Since this
    counter is based on a tradeoff between precision and memory consumption,
    the :option:`-e`/:option:`--error-rate` can be used to control how much
    memory will be used. In practice the memory footprint is small even
    at low error rates (< 0.01).

    :option:`-k`/:option:`--ksize` should be set to the desired k-mer size.

    Informational output is sent to STDERR, but a report file can be generated
    with :option:`-R`/:option:`--report`.

    :option:`--stream-records` will write the sequences taken in to STDOUT.
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

        unique-kmers.py --stream-records -k 17 tests/test-data/test-reads.fa | \\
        normalize-by-median.py -k 17 -o normalized /dev/stdin

    Example::

        unique-kmers.py -R unique_count -k 30 \\
        tests/test-data/test-abund-read-paired.fa"""  # noqa
    parser = KhmerArgumentParser(
        description=descr, epilog=textwrap.dedent(epilog),
        citations=['SeqAn', 'hll'])

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)

    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')

    parser.add_argument('-k', '--ksize', type=int, default=env_ksize,
                        help='k-mer size to use')

    parser.add_argument('-e', '--error-rate', type=float, default=0.01,
                        help='Acceptable error rate')

    parser.add_argument('-R', '--report',
                        metavar='filename', type=argparse.FileType('w'),
                        help='generate informational report and write to'
                        ' filename')

    parser.add_argument('-S', '--stream-records', default=False,
                        action='store_true',
                        help='write input sequences to STDOUT')

    parser.add_argument('--diagnostics', default=False, action='store_true',
                        help='print out recommended tablesize arguments and '
                             'restrictions')

    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename(s).', nargs='+')

    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()

    total_hll = khmer.HLLCounter(args.error_rate, args.ksize)

    report_fp = args.report
    input_filename = None
    for _, input_filename in enumerate(args.input_filenames):
        hllcpp = khmer.HLLCounter(args.error_rate, args.ksize)
        hllcpp.consume_seqfile(input_filename,
                               stream_records=args.stream_records)

        cardinality = hllcpp.estimate_cardinality()
        print('Estimated number of unique {0}-mers in {1}: {2}'.format(
            args.ksize, input_filename, cardinality), file=sys.stderr)

        if report_fp:
            print(cardinality, args.ksize, '(total)', file=report_fp)
            report_fp.flush()
        total_hll.merge(hllcpp)

    cardinality = total_hll.estimate_cardinality()
    print('Total estimated number of unique {0}-mers: {1}'.format(
        args.ksize, cardinality), file=sys.stderr)

    to_print = graphsize_args_report(cardinality, args.error_rate)
    if args.diagnostics:
        print(to_print, file=sys.stderr)

    if report_fp:
        print(cardinality, args.ksize, 'total', file=report_fp)
        print(to_print, file=report_fp)
        report_fp.flush()


if __name__ == "__main__":
    main()
