#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring

"""
Convert FASTQ files to FASTA format.

% python scripts/fastq-to-fasta.py [ -n -o ] <fastq_name>

Use '-h' for parameter help.
"""
import sys
import argparse
import screed


def get_parser():
    parser = argparse.ArgumentParser(
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_sequence', help='The name of the input'
                        ' FASTQ sequence file.')
    parser.add_argument('-o', '--output', help='The name of the output'
                        ' FASTA sequence file.')
    parser.add_argument('-n', '--n_keep', default=False, action='store_true',
                        help='Option to drop reads containing \'N\'s in ' +
                        'input_sequence file.')
    return parser


def main():
    args = get_parser().parse_args()
    print >> sys.stderr, ('fastq from ', args.input_sequence)

    if args.output:
        write_out = open(args.output, 'w')

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
        if args.output:
            write_out.write('>' + name + '\n')
            write_out.write(sequence + '\n')
        else:
            print '>' + name
            print sequence

    print >> sys.stderr, '\n' + 'lines from ' + args.input_sequence

    if not args.n_keep:
        print >> sys.stderr, str(n_count) + ' lines dropped.'

    else:
        print >> sys.stderr, 'No lines dropped from file.'

if __name__ == '__main__':
    main()
