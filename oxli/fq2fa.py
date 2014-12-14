#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring

"""
Convert FASTQ files to FASTA format.
"""
import sys
import argparse
import screed


def add_args(parser):
    parser.add_argument('input_sequence', help='The name of the input'
                              ' FASTQ sequence file.')
    parser.add_argument('-o', '--output', metavar="filename",
                              help='The name of the output'
                              ' FASTA sequence file.',
                              type=argparse.FileType('w'),
                              default=sys.stdout)
    parser.add_argument('-n', '--n_keep', default=False,
                              action='store_true', help='Option to drop reads containing \'N\'s \
            in input_sequence file.')


def do_fastq_to_fasta(args):
    #args = get_parser().parse_args()
    print >> sys.stderr, ('fastq from ', args.input_sequence)

    n_count = 0
    for n, record in enumerate(screed.open(args.input_sequence)):
        if n % 10000 == 0:
            print>>sys.stderr, '...', n

        sequence = record['sequence']
        name = record['name']

        if 'N' in sequence:
            if not args.n_keep:
                n_count += 1
                continue

        args.output.write('>' + name + '\n')
        args.output.write(sequence + '\n')

    print >> sys.stderr, '\n' + 'lines from ' + args.input_sequence

    if not args.n_keep:
        print >> sys.stderr, str(n_count) + ' lines dropped.'

    else:
        print >> sys.stderr, 'No lines dropped from file.'

    print >> sys.stderr, 'Wrote output to', args.output

