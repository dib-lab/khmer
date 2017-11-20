#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
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
Subsample sequences from multiple files.

Take a list of files containing sequences, and subsample 100,000 sequences (-N)
uniformly, using reservoir sampling.  Stop after first 100m sequences (-M).
By default take one subsample, but take -S samples if specified.

% scripts/sample-reads-randomly.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""

import argparse
import os.path
import random
import textwrap
import sys

from khmer import __version__
from khmer import ReadParser
from khmer.kfile import (check_input_files, add_output_compression_type,
                         get_file_writer)
from khmer.khmer_args import sanitize_help, KhmerArgumentParser
from khmer.utils import write_record, broken_paired_reader

DEFAULT_NUM_READS = int(1e5)
DEFAULT_MAX_READS = int(1e8)
DEBUG = True


def get_parser():
    epilog = """\
    Take a list of files containing sequences, and subsample 100,000
    sequences (:option:`-N`/:option:`--num_reads`) uniformly, using
    reservoir sampling.  Stop after first 100m sequences
    (:option:`-M`/:option:`--max_reads`). By default take one subsample,
    but take :option:`-S`/:option:`--samples` samples if specified.

    The output is placed in :option:`-o`/:option:`--output` <file>
    (for a single sample) or in ``<file>.subset.0`` to ``<file>.subset.S-1``
    (for more than one sample).

    This script uses the `reservoir sampling
    <http://en.wikipedia.org/wiki/Reservoir_sampling>`__ algorithm.
    """

    parser = KhmerArgumentParser(
        description="Uniformly subsample sequences from a collection of files",
        epilog=textwrap.dedent(epilog))

    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-N', '--num_reads', type=int, dest='num_reads',
                        default=DEFAULT_NUM_READS, help='samples the '
                        'number of sequences or pairs specified with -N')
    parser.add_argument('-M', '--max_reads', type=int, dest='max_reads',
                        default=DEFAULT_MAX_READS)
    parser.add_argument('-S', '--samples', type=int, dest='num_samples',
                        default=1)
    parser.add_argument('-R', '--random-seed', type=int, dest='random_seed',
                        help='Provide a random seed for the generator')
    parser.add_argument('--force_single', default=False, action='store_true',
                        help='Ignore read pair information if present')
    parser.add_argument('-o', '--output', dest='output_file',
                        type=argparse.FileType('wb'),
                        metavar="filename", default=None)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exits')
    add_output_compression_type(parser)
    return parser


def main():
    parser = get_parser()
    parser.epilog = parser.epilog.replace(
        "`reservoir sampling\n"
        "<http://en.wikipedia.org/wiki/Reservoir_sampling>`__ algorithm.",
        "reservoir sampling algorithm. "
        "http://en.wikipedia.org/wiki/Reservoir_sampling")
    args = sanitize_help(parser).parse_args()

    for name in args.filenames:
        check_input_files(name, args.force)

    # seed the random number generator?
    if args.random_seed:
        random.seed(args.random_seed)

    # bound n_samples
    num_samples = max(args.num_samples, 1)

    #
    # Figure out what the output filename is going to be

    if args.output_file:
        output_filename = args.output_file.name
        if num_samples > 1:
            sys.stderr.write(
                "Error: cannot specify -o with more than one sample.")
            if not args.force:
                print("NOTE: This can be overridden using the --force"
                      " argument", file=sys.stderr)
                sys.exit(1)
    else:
        filename = args.filenames[0]
        if filename in ('/dev/stdin', '-'):
            print("Accepting input from stdin; output filename must "
                  "be provided with '-o'.", file=sys.stderr)
            sys.exit(1)
        output_filename = os.path.basename(filename) + '.subset'

    filename = args.filenames[0]
    if filename in ('/dev/stdin', '-'):
        # seqan only treats '-' as "read from stdin"
        filename = '-'

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
    for _ in range(num_samples):
        reads.append([])

    # read through all the sequences and load/resample the reservoir
    for filename in args.filenames:
        print('opening', filename, 'for reading', file=sys.stderr)

        for count, (_, _, rcrd1, rcrd2) in enumerate(broken_paired_reader(
                ReadParser(filename), force_single=args.force_single)):
            if count % 10000 == 0:
                print('...', count, 'reads scanned', file=sys.stderr)
                if count >= args.max_reads:
                    print('reached upper limit of %d reads' %
                          args.max_reads, '(see -M); exiting', file=sys.stderr)
                    break

            # collect first N reads
            if count < args.num_reads:
                for sample in range(num_samples):
                    reads[sample].append((rcrd1, rcrd2))
            else:
                for sample in range(num_samples):
                    assert len(reads[sample]) <= count

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

        output_file = args.output_file
        if not output_file:
            output_file = open(output_filename, 'wb')

        output_file = get_file_writer(output_file, args.gzip, args.bzip)

        for records in reads[0]:
            write_record(records[0], output_file)
            if records[1] is not None:
                write_record(records[1], output_file)
    else:
        for n in range(num_samples):
            n_filename = output_filename + '.%d' % n
            print('Writing %d sequences to %s' %
                  (len(reads[n]), n_filename), file=sys.stderr)
            output_file = get_file_writer(open(n_filename, 'wb'), args.gzip,
                                          args.bzip)
            for records in reads[n]:
                write_record(records[0], output_file)
                if records[1] is not None:
                    write_record(records[1], output_file)


if __name__ == '__main__':
    main()
