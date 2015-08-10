#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Produce the median k-mer abundance distribution for reads in the given file.
"""
from __future__ import division
from __future__ import print_function

import sys
import argparse
import os
import screed
import khmer
import csv
from khmer.kfile import check_input_files
from khmer.khmer_args import info

# TODO: rename to median-dist.py => scripts/.

def get_parser():
    parser = argparse.ArgumentParser(
        description="Output the distribution of median k-mer abundances per "
        "read, which is a reference-free proxy for read coverage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_countgraph', help='The name of the input '
                        'k-mer countgraph file, e.g. produced by '
                        'load-into-counting.py')
    parser.add_argument('input_sequence_filename', help='The name of the '
                        'FAST[AQ] sequence file.')
    parser.add_argument('output_histogram_filename', help='@@')

    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite existing output_histogram_filename')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count abundances past 255')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Continue even if specified input files '
                        'do not exist or are empty.')

    return parser


def main():
    info('abundance-dist.py', ['counting'])
    args = get_parser().parse_args()

    infiles = [args.input_countgraph,
               args.input_sequence_filename]
    for infile in infiles:
        check_input_files(infile, False)

    print('Counting graph from', args.input_countgraph,
          file=sys.stderr)
    countgraph = khmer.load_countgraph(args.input_countgraph)

    if not countgraph.get_use_bigcount() and args.bigcount:
        print("WARNING: The loaded graph has bigcount DISABLED while bigcount"
              " reporting is ENABLED--counts higher than 255 will not be "
              "reported.",
              file=sys.stderr)

    countgraph.set_use_bigcount(args.bigcount)

    # @CTB
    kmer_size = countgraph.ksize()
    hashsizes = countgraph.hashsizes()
    tracking = khmer._Nodegraph(  # pylint: disable=protected-access
        kmer_size, hashsizes)

    print('K:', kmer_size, file=sys.stderr)
    print('outputting to', args.output_histogram_filename, file=sys.stderr)

    if args.output_histogram_filename in ('-', '/dev/stdout'):
        pass
    elif os.path.exists(args.output_histogram_filename):
        if not args.squash_output:
            print('ERROR: %s exists; not squashing.' %
                  args.output_histogram_filename,
                  file=sys.stderr)
            sys.exit(1)

        print('** squashing existing file %s' %
              args.output_histogram_filename, file=sys.stderr)

    # prepare histogram
    hist = {}
    for i in range(65536):
        hist[i] = 0

    # walk through the reads and count.
    for n, record in enumerate(screed.open(args.input_sequence_filename)):
        if n > 0 and n % 100000 == 0:
            print('...', n)

        seq = record.sequence.replace('N', 'A')

        try:
            med, _, _ = countgraph.get_median_count(seq)
        except ValueError:
            continue

        hist[med] = hist[med] + 1

    # ok, we got all the median abundances in 'hist'; summarize & output
    # the distribution.
    total = sum(hist.values())
    abundances = list(hist.items())
    abundances.sort()

    if 0 == total:
        print("ERROR: abundance distribution is uniformly zero; "
              "nothing to report.", file=sys.stderr)
        print("\tPlease verify that the input files are valid.",
              file=sys.stderr)
        sys.exit(1)

    # pick output.
    if args.output_histogram_filename in ('-', '/dev/stdout'):
        countgraph_fp = sys.stdout
    else:
        countgraph_fp = open(args.output_histogram_filename, 'w')
    countgraph_fp_csv = csv.writer(countgraph_fp)
    # write headers:
    countgraph_fp_csv.writerow(['abundance', 'count', 'cumulative',
                                'cumulative_fraction'])

    # produce output.
    sofar = 0
    for n, i in abundances:
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)
        countgraph_fp_csv.writerow([n, i, sofar, round(frac, 3)])

        if sofar == total:
            break


if __name__ == '__main__':
    main()
