#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
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
# pylint: disable=invalid-name,missing-docstring
"""
Produce the k-mer abundance distribution for the given file.

% python scripts/abundance-dist-single.py <data> <histout>

The script does not load a prebuilt k-mer countgraph.

Use '-h' for parameter help.
"""
from __future__ import print_function
import os
import sys
import csv
import khmer
import threading
import textwrap
from khmer import khmer_args
from khmer.khmer_args import (build_counting_args, add_threading_args,
                              report_on_config, info, calculate_graphsize,
                              sanitize_help)
from khmer.kfile import (check_input_files, check_space_for_graph)


def get_parser():
    epilog = '''\
    Note that with :option:`-b`/:option:`--no-bigcount` this script is constant
    memory; in exchange, k-mer counts will stop at 255. The memory usage of
    this script with :option:`-b` will be about 1.15x the product of the
    :option:`-x` and :option:`-N` numbers.

    To count k-mers in multiple files use :program:`load_into_counting.py` and
    :program:`abundance_dist.py`.

    Example::

        abundance-dist-single.py -x 1e7 -N 2 -k 17 \\
                tests/test-data/test-abund-read-2.fa test-dist
    '''
    parser = build_counting_args(
        descr="Calculate the abundance distribution of k-mers from a "
        "single sequence file.", epilog=textwrap.dedent(epilog))
    add_threading_args(parser)

    parser.add_argument('input_sequence_filename', help='The name of the input'
                        ' FAST[AQ] sequence file.')
    parser.add_argument('output_histogram_filename', help='The name of the '
                        'output histogram file. The columns are: (1) k-mer '
                        'abundance, (2) k-mer count, (3) cumulative count, '
                        '(4) fraction of total distinct k-mers.')
    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output zero-count bins')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--savegraph', default='', metavar="filename",
                        help="Save the k-mer countgraph to the specified "
                        "filename.")
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
    return parser


def main():  # pylint: disable=too-many-locals,too-many-branches
    args = sanitize_help(get_parser()).parse_args()
    if not args.quiet:
        info('abundance-dist-single.py', ['counting', 'SeqAn'])
    report_on_config(args)

    check_input_files(args.input_sequence_filename, args.force)
    if args.savegraph:
        graphsize = calculate_graphsize(args, 'countgraph')
        check_space_for_graph(args.savegraph, graphsize, args.force)
    if (not args.squash_output and
            os.path.exists(args.output_histogram_filename)):
        print('ERROR: %s exists; not squashing.' %
              args.output_histogram_filename, file=sys.stderr)
        sys.exit(1)
    else:
        hist_fp = open(args.output_histogram_filename, 'w')
        hist_fp_csv = csv.writer(hist_fp)
        # write headers:
        hist_fp_csv.writerow(['abundance', 'count', 'cumulative',
                              'cumulative_fraction'])

    print('making countgraph', file=sys.stderr)
    countgraph = khmer_args.create_countgraph(args, multiplier=1.1)
    countgraph.set_use_bigcount(args.bigcount)

    print('building k-mer tracking graph', file=sys.stderr)
    tracking = khmer_args.create_nodegraph(args, multiplier=1.1)

    print('kmer_size:', countgraph.ksize(), file=sys.stderr)
    print('k-mer countgraph sizes:',
          countgraph.hashsizes(), file=sys.stderr)
    print('outputting to', args.output_histogram_filename, file=sys.stderr)

    # start loading
    rparser = khmer.ReadParser(args.input_sequence_filename)
    threads = []
    print('consuming input, round 1 --',
          args.input_sequence_filename, file=sys.stderr)
    for _ in range(args.threads):
        thread = \
            threading.Thread(
                target=countgraph.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    print('Total number of unique k-mers: {0}'.format(
        countgraph.n_unique_kmers()), file=sys.stderr)

    abundance_lists = []

    def __do_abundance_dist__(read_parser):
        abundances = countgraph.abundance_distribution_with_reads_parser(
            read_parser, tracking)
        abundance_lists.append(abundances)

    print('preparing hist from %s...' %
          args.input_sequence_filename, file=sys.stderr)
    rparser = khmer.ReadParser(args.input_sequence_filename)
    threads = []
    print('consuming input, round 2 --',
          args.input_sequence_filename, file=sys.stderr)
    for _ in range(args.threads):
        thread = \
            threading.Thread(
                target=__do_abundance_dist__,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    assert len(abundance_lists) == args.threads, len(abundance_lists)
    abundance = {}
    for abundance_list in abundance_lists:
        for i, count in enumerate(abundance_list):
            abundance[i] = abundance.get(i, 0) + count

    total = sum(abundance.values())

    if 0 == total:
        print("ERROR: abundance distribution is uniformly zero; "
              "nothing to report.", file=sys.stderr)
        print(
            "\tPlease verify that the input files are valid.", file=sys.stderr)
        sys.exit(1)

    sofar = 0
    for _, i in sorted(abundance.items()):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        hist_fp_csv.writerow([_, i, sofar, round(frac, 3)])

        if sofar == total:
            break

    if args.savegraph:
        print('Saving k-mer countgraph ', args.savegraph, file=sys.stderr)
        print('...saving to', args.savegraph, file=sys.stderr)
        countgraph.save(args.savegraph)

    print('wrote to: ' + args.output_histogram_filename, file=sys.stderr)

if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
