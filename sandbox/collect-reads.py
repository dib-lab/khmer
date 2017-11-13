#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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
# pylint: disable=missing-docstring,invalid-name
"""
Build a counting Bloom filter from the given sequences, save in <htname>.
Stop collecting reads when the average coverage gets above -C (default 50).
Place reads into -o output_file.

% collect-reads.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import textwrap
import khmer
from khmer import khmer_args
from khmer.khmer_args import (build_counting_args, report_on_config, info,
                              calculate_graphsize, sanitize_help)
from khmer.kfile import check_input_files, check_space
from khmer.kfile import check_space_for_graph
import argparse
import screed


def output_single(read):
    if hasattr(read, 'quality'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.quality)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def get_parser():
    epilog = """
    The memory usage of this script with :option:`-b` will be about
    1.15x the product of the :option:`-x` and :option:`-N` numbers.

    Example::

        collect-reads.py -k 20 -x 5e7 out.ct data/100k-filtered.fa
    """

    parser = build_counting_args("Collect reads until a given avg coverage.",
                                 epilog=textwrap.dedent(epilog))
    parser.add_argument('output_countgraph_filename', help="The name of the"
                        " file to write the k-mer countgraph to.")
    parser.add_argument('input_sequence_filename', nargs='+',
                        help="The names of one or more FAST[AQ] input "
                        "sequence files.")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    parser.add_argument('-C', '--coverage', type=int, default=50,
                        help='Collect reads until this coverage, then exit.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Write collect reads into this file.')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    return parser


def main():

    info('collect-reads.py', ['counting'])
    args = sanitize_help(get_parser()).parse_args()
    report_on_config(args)

    base = args.output_countgraph_filename
    filenames = args.input_sequence_filename

    for name in args.input_sequence_filename:
        check_input_files(name, False)

    check_space(args.input_sequence_filename, False)
    tablesize = calculate_graphsize(args, 'countgraph')
    check_space_for_graph(args.output_countgraph_filename, tablesize,
                              False)

    print('Saving k-mer countgraph to %s' % base)
    print('Loading sequences from %s' % repr(filenames))
    if args.output:
        print('Outputting sequences to', args.output)

    print('making countgraph', file=sys.stderr)
    htable = khmer_args.create_countgraph(args)
    htable.set_use_bigcount(args.bigcount)

    total_coverage = 0.
    n = 0

    for index, filename in enumerate(filenames):
        for record in screed.open(filename):
            seq = record.sequence.upper()
            if 'N' in seq:
                seq = seq.replace('N', 'A')

            try:
                med, _, _ = htable.get_median_count(seq)
            except ValueError:
                continue

            total_coverage += med
            n += 1

            if total_coverage / float(n) > args.coverage:
                print('reached target average coverage:', \
                      total_coverage / float(n))
                break

            htable.consume(seq)
            if args.output:
                args.output.write(output_single(record))

            if n % 100000 == 0:
                print('...', index, filename, n, total_coverage / float(n))

        if total_coverage / float(n) > args.coverage:
            break

    print('Collected %d reads' % (n,))

    if args.report_total_kmers:
        print('Total number of k-mers: {0}'.format(
            htable.n_occupied()), file=sys.stderr)

    print('saving', base)
    htable.save(base)

    info_fp = open(base + '.info', 'w')
    info_fp.write('through end: %s\n' % filenames[-1])

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(htable, False,
                                             max_false_pos=.2)
    print('fp rate estimated to be %1.3f' % fp_rate)
    print('fp rate estimated to be %1.3f' % fp_rate, file=info_fp)

    print('DONE.')

if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
