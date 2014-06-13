#! /usr/bin/env python2
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Take a list of files containing sequences, and subsample 100,000 sequences (-N)
uniformly, using reservoir sampling.  Stop after first 100m sequences (-M).
By default take one subsample, but take -S samples if specified.

% scripts/sample-reads-randomly.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""

import argparse
import screed
import os.path
import random
import textwrap
import sys

import khmer
from khmer.file import check_file_status, check_space
from khmer.khmer_args import info

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
    parser.add_argument('-o', '--output', dest='output_file',
                        metavar='output_file',
                        type=argparse.FileType('w'), default=None)
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def output_single(read):
    if hasattr(read, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.accuracy)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def main():
    info('sample-reads-randomly.py')
    args = get_parser().parse_args()

    for _ in args.filenames:
        check_file_status(_)

    check_space(args.filenames)

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
            sys.exit(-1)
        output_filename = output_file.name
    else:
        filename = args.filenames[0]
        output_filename = os.path.basename(filename) + '.subset'

    if num_samples == 1:
        print 'Subsampling %d reads using reservoir sampling.' % args.num_reads
        print 'Subsampled reads will be placed in %s' % output_filename
        print ''
    else:  # > 1
        print 'Subsampling %d reads, %d times, using reservoir sampling.' % \
            (args.num_reads, num_samples)
        print 'Subsampled reads will be placed in %s.N' % output_filename
        print ''

    reads = []
    for n in range(num_samples):
        reads.append([])

    total = 0

    # read through all the sequences and load/resample the reservoir
    for filename in args.filenames:
        print 'opening', filename, 'for reading'
        for record in screed.open(filename):
            total += 1

            if total % 10000 == 0:
                print '...', total, 'reads scanned'
                if total >= args.max_reads:
                    print 'reached upper limit of %d reads (see -M); exiting' \
                        % args.max_reads
                    break

            # collect first N reads
            if total <= args.num_reads:
                for n in range(num_samples):
                    reads[n].append(record)
            else:
                # use reservoir sampling to replace reads at random
                # see http://en.wikipedia.org/wiki/Reservoir_sampling

                for n in range(num_samples):
                    guess = random.randint(1, total)
                    if guess <= args.num_reads:
                        reads[n][guess - 1] = record

    # output all the subsampled reads:
    if len(reads) == 1:
        print 'Writing %d sequences to %s' % (len(reads[0]), output_filename)
        if not output_file:
            output_file = open(output_filename, 'w')

        for record in reads[0]:
            output_file.write(output_single(record))
    else:
        for n in range(num_samples):
            n_filename = output_filename + '.%d' % n
            print 'Writing %d sequences to %s' % (len(reads[n]), n_filename)
            output_file = open(n_filename, 'w')
            for record in reads[n]:
                output_file.write(output_single(record))

if __name__ == '__main__':
    main()
