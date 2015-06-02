#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Produce the k-mer abundance distribution for the given file.

% python scripts/abundance-dist-single.py <data> <histout>

The script does not load a prebuilt k-mer counting table.

Use '-h' for parameter help.
"""
from __future__ import print_function
import os
import sys
import csv
import khmer
import threading
import textwrap
from khmer.khmer_args import (build_counting_args, add_threading_args,
                              report_on_config, info)
from khmer.kfile import (check_input_files, check_space,
                         check_space_for_hashtable)


def get_parser():
    epilog = '''
    Note that with :option:`-b` this script is constant memory; in exchange,
    k-mer counts will stop at 255. The memory usage of this script with
    :option:`-b` will be about 1.15x the product of the :option:`-x` and
    :option:`-N` numbers.

    To count k-mers in multiple files use :program:`load_into_counting.py` and
    :program:`abundance_dist.py`.
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
                        help='Do not output 0-count bins')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--csv', default=False, action='store_true',
                        help='Use the CSV format for the histogram. '
                        'Includes column headers.')
    parser.add_argument('--savetable', default='', metavar="filename",
                        help="Save the k-mer counting table to the specified "
                        "filename.")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():  # pylint: disable=too-many-locals,too-many-branches
    info('abundance-dist-single.py', ['counting', 'SeqAn'])
    args = get_parser().parse_args()
    report_on_config(args)

    check_input_files(args.input_sequence_filename, args.force)
    check_space([args.input_sequence_filename], args.force)
    if args.savetable:
        check_space_for_hashtable(args.n_tables * args.min_tablesize,
                                  args.force)

    if (not args.squash_output and
            os.path.exists(args.output_histogram_filename)):
        print('ERROR: %s exists; not squashing.' %
              args.output_histogram_filename, file=sys.stderr)
        sys.exit(1)
    else:
        hist_fp = open(args.output_histogram_filename, 'w')
        if args.csv:
            hist_fp_csv = csv.writer(hist_fp)
            # write headers:
            hist_fp_csv.writerow(['abundance', 'count', 'cumulative',
                                  'cumulative_fraction'])

    print('making k-mer counting table', file=sys.stderr)
    counting_hash = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                            args.n_tables)
    counting_hash.set_use_bigcount(args.bigcount)

    print('building k-mer tracking table', file=sys.stderr)
    tracking = khmer.new_hashbits(counting_hash.ksize(), args.min_tablesize,
                                  args.n_tables)

    print('kmer_size:', counting_hash.ksize(), file=sys.stderr)
    print('k-mer counting table sizes:',
          counting_hash.hashsizes(), file=sys.stderr)
    print('outputting to', args.output_histogram_filename, file=sys.stderr)

    # start loading
    rparser = khmer.ReadParser(args.input_sequence_filename)
    threads = []
    print('consuming input, round 1 --',
          args.input_sequence_filename, file=sys.stderr)
    for _ in range(args.threads):
        thread = \
            threading.Thread(
                target=counting_hash.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    if args.report_total_kmers:
        print('Total number of unique k-mers: {0}'.format(
            counting_hash.n_unique_kmers()), file=sys.stderr)

    abundance_lists = []

    def __do_abundance_dist__(read_parser):
        abundances = counting_hash.abundance_distribution_with_reads_parser(
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

        if args.csv:
            hist_fp_csv.writerow([_, i, sofar, round(frac, 3)])
        else:
            print(_, i, sofar, round(frac, 3), file=hist_fp)

        if sofar == total:
            break

    if args.savetable:
        print('Saving k-mer counting table ', args.savetable, file=sys.stderr)
        print('...saving to', args.savetable, file=sys.stderr)
        counting_hash.save(args.savetable)

    print('wrote to: ' + args.output_histogram_filename, file=sys.stderr)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
