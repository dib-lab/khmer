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
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_sequence', help='The name of the input'
                        ' FASTQ sequence file.')
    parser.add_argument('-o', '--output', help='The name of the output'
                        ' FASTA sequence file.')
    return parser


def main():
    min_length = int(sys.argv[1])
    if args.output:
        write_out = open(args.output, 'w')

    for filename in sys.argv[2:]:
        for record in screed.open(filename):
            if len(record['sequence']) >= min_length:
                #if fasta

                    print >> sys.stderr, ( '>%s\n%s' % (record['name'], 
                        record['sequence'],) )
                #elif fastq
                    if hasattr(record, 'accuracy'):
                        outfp.write(
                            '@{name}\n{seq}\n'
                            '+\n{acc}\n'.format(name=record.name,
                                                seq=record.sequence,
                                                acc=record.accuracy))
                    else:
                        outfp.write(
                            '>{name}\n{seq}\n'.format(name=record.name,
                                                        seq=record.sequence))

if __name__ == '__main__':
    main()
