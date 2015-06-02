#! /usr/bin/env python
#
# This script is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Subsample sequences from multiple files.

Take a list of files containing sequences, and subsample 100,000 sequences (-N)
uniformly, using reservoir sampling.  Stop after first 100m sequences (-M).
By default take one subsample, but take -S samples if specified.

% scripts/sample-reads-randomly.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""
from __future__ import print_function

import argparse
import screed
import os.path
import random
import textwrap
import sys

import khmer
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info
from khmer.utils import write_record, broken_paired_reader

DEFAULT_NUM_READS = int(1e5)
DEFAULT_MAX_READS = int(1e8)
DEBUG = True


def get_parser():
    epilog = ("""

    Take a list of files containing sequences, and subsample 100,000
    sequences (:option:`-N`/:option:`--num_reads`) uniformly, using
    reservoir sampling.  Stop after first 100m sequences
    (:option:`-M`/:option:`--max_reads`). By default take one subsample,
    but take :option:`-S`/:option:`--samples` samples if specified.

    The output is placed in :option:`-o`/:option:`--output` <file>
    (for a single sample) or in <file>.subset.0 to <file>.subset.S-1
    (for more than one sample).

    This script uses the `reservoir sampling
    <http://en.wikipedia.org/wiki/Reservoir_sampling>`__ algorithm.
    """)   # noqa

    parser = argparse.ArgumentParser(
        description="Uniformly subsample sequences from a collection of files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=textwrap.dedent(epilog))

    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-N', '--num_reads', type=int, dest='num_reads',
                        default=DEFAULT_NUM_READS)
    parser.add_argument('-M', '--max_reads', type=int, dest='max_reads',
                        default=DEFAULT_MAX_READS)
    parser.add_argument('-S', '--samples', type=int, dest='num_samples',
                        default=1)
    parser.add_argument('-R', '--random-seed', type=int, dest='random_seed')
    parser.add_argument('--force_single', default=False, action='store_true',
                        help='Ignore read pair information if present')
    parser.add_argument('-o', '--output', dest='output_file',
                        metavar='output_file',
                        type=argparse.FileType('w'), default=None)
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exits')
    return parser


def main():
    info('sample-reads-randomly.py')
    args = get_parser().parse_args()

    for _ in args.filenames:
        check_input_files(_, args.force)

    check_space(args.filenames, args.force)

    # seed the random number generator?
    if args.random_seed:
        random.seed(args.random_seed)

    # bound n_samples
    num_samples = max(args.num_samples, 1)

    #
    # Figure out what the output filename is going to be
    #

    output_file = args.output_file
    if output_file:
        if num_samples > 1:
            sys.stderr.write(
                "Error: cannot specify -o with more than one sample.")
            if not args.force:
                sys.exit(1)
        output_filename = output_file.name
    else:
        filename = args.filenames[0]
        output_filename = os.path.basename(filename) + '.subset'

    if num_samples == 1:
        print('Subsampling %d reads using reservoir sampling.' %
              args.num_reads, file=sys.stderr)
        print('Subsampled reads will be placed in %s' %
              output_filename, file=sys.stderr)
        print('', file=sys.stderr)
    else:  # > 1
        print('Subsampling %d reads, %d times,'
              % (args.num_reads, num_samples), ' using reservoir sampling.',
              file=sys.stderr)
        print('Subsampled reads will be placed in %s.N'
              % output_filename, file=sys.stderr)
        print('', file=sys.stderr)

    reads = []
    for n in range(num_samples):
        reads.append([])

    # read through all the sequences and load/resample the reservoir
    for filename in args.filenames:
        print('opening', filename, 'for reading', file=sys.stderr)
        screed_iter = screed.open(filename, parse_description=False)

        for count, (_, ispair, rcrd1, rcrd2) in enumerate(broken_paired_reader(
                screed_iter,
                force_single=args.force_single)):
            if count % 10000 == 0:
                print('...', count, 'reads scanned', file=sys.stderr)
                if count >= args.max_reads:
                    print('reached upper limit of %d reads' %
                          args.max_reads, '(see -M); exiting', file=sys.stderr)
                    break

            # collect first N reads
            if count < args.num_reads:
                for n in range(num_samples):
                    reads[n].append((rcrd1, rcrd2))
            else:
                assert len(reads[n]) <= count

                # use reservoir sampling to replace reads at random
                # see http://en.wikipedia.org/wiki/Reservoir_sampling

                for n in range(num_samples):
                    guess = random.randint(1, count)
                    if guess <= args.num_reads:
                        reads[n][guess - 1] = (rcrd1, rcrd2)

    # output all the subsampled reads:
    if len(reads) == 1:
        print('Writing %d sequences to %s' %
              (len(reads[0]), output_filename), file=sys.stderr)
        if not output_file:
            output_file = open(output_filename, 'w')

        for records in reads[0]:
            write_record(records[0], output_file)
            if records[1] is not None:
                write_record(records[1], output_file)
    else:
        for n in range(num_samples):
            n_filename = output_filename + '.%d' % n
            print('Writing %d sequences to %s' %
                  (len(reads[n]), n_filename), file=sys.stderr)
            output_file = open(n_filename, 'w')
            for records in reads[n]:
                write_record(records[0], output_file)
                if records[1] is not None:
                    write_record(records[1], output_file)

if __name__ == '__main__':
    main()
