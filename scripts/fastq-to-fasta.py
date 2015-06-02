#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring

"""
Convert FASTQ files to FASTA format.

% python scripts/fastq-to-fasta.py [ -n -o ] <fastq_name>

Use '-h' for parameter help.
"""
from __future__ import print_function
import sys
import argparse
import screed


def get_parser():
    parser = argparse.ArgumentParser(
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_sequence', help='The name of the input'
                        ' FASTQ sequence file.')
    parser.add_argument('-o', '--output', metavar="filename",
                        help='The name of the output'
                        ' FASTA sequence file.',
                        type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-n', '--n_keep', default=False, action='store_true',
                        help='Option to drop reads containing \'N\'s in ' +
                        'input_sequence file.')
    return parser


def main():
    args = get_parser().parse_args()
    print(('fastq from ', args.input_sequence), file=sys.stderr)

    n_count = 0
    for n, record in enumerate(screed.open(args.input_sequence,
                                           parse_description=False)):
        if n % 10000 == 0:
            print('...', n, file=sys.stderr)

        sequence = record['sequence']
        name = record['name']

        if 'N' in sequence:
            if not args.n_keep:
                n_count += 1
                continue

        args.output.write('>' + name + '\n')
        args.output.write(sequence + '\n')

    print('\n' + 'lines from ' + args.input_sequence, file=sys.stderr)

    if not args.n_keep:
        print(str(n_count) + ' lines dropped.', file=sys.stderr)

    else:
        print('No lines dropped from file.', file=sys.stderr)

    print('Wrote output to', args.output, file=sys.stderr)

if __name__ == '__main__':
    main()
