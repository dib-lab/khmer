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
import argparse
import screed
import sys


def get_parser():
    parser = argparse.ArgumentParser(
        description='Extracts FASTQ or FASTA sequences longer than argument 
        specified length.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_filenames', help='Input FAST[AQ]'
                        ' sequence filename.', nargs='+')
    parser.add_argument('-o', '--output', help='The name of the output'
                        ' sequence file.')
    parser.add_argument('-l', '--length', help='The minimum length of the'
                        ' sequence file.')
    return parser


def main():
    args = get_parser().parse_args()
    if args.output:
        write_out = open(args.output, 'w')

    for file in args.input_filenames
        for record in screed.open(file):
            if len(record['sequence']) >= args.length:
                if ('+' in file):
                    print >> sys.stderr, ( '>%s\n%s' % (record['name'], 
                        record['sequence'],) )
                else:
                    if hasattr(record, 'accuracy'):
                        write_out.write(
                            '@{name}\n{seq}\n'
                            '+\n{acc}\n'.format(name=record.name,
                                                seq=record.sequence,
                                                acc=record.accuracy))
                    else:
                        write_out.write(
                            '>{name}\n{seq}\n'.format(name=record.name,
                                                        seq=record.sequence))

if __name__ == '__main__':
    main()
