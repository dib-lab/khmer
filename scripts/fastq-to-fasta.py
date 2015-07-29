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
from __future__ import print_function, unicode_literals
import sys
import argparse
import screed
from khmer.kfile import (add_output_compression_type, get_file_writer,
                         is_block)
from khmer.utils import write_record


def get_parser():
    parser = argparse.ArgumentParser(
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_sequence', help='The name of the input'
                        ' FASTQ sequence file.')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=argparse.FileType('wb'),
                        help='The name of the output'
                        ' FASTA sequence file.',
                        default=sys.stdout)
    parser.add_argument('-n', '--n_keep', default=False, action='store_true',
                        help='Option to keep reads containing \'N\'s in '
                             'input_sequence file. Default is to drop reads')
    add_output_compression_type(parser)
    return parser


def main():
    args = get_parser().parse_args()

    print(('fastq from ', args.input_sequence), file=sys.stderr)
    output_is_block = is_block(args.output)
    outfp = get_file_writer(args.output, args.gzip, args.bzip)
    n_count = 0
    for n, record in enumerate(screed.open(args.input_sequence)):
        if n % 10000 == 0:
            print('...', n, file=sys.stderr)

        sequence = record['sequence']
        name = record['name']

        if 'N' in sequence:
            if not args.n_keep:
                n_count += 1
                continue

        del record['quality']
        write_record(record, outfp)

    print('\n' + 'lines from ' + args.input_sequence, file=sys.stderr)

    if not args.n_keep:
        print(str(n_count) + ' lines dropped.', file=sys.stderr)

    else:
        print('No lines dropped from file.', file=sys.stderr)

    if output_is_block:
        print('Wrote output to block device', file=sys.stderr)
    else:
        print('Wrote output to', args.output.name, file=sys.stderr)

if __name__ == '__main__':
    main()
