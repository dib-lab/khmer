#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring

"""
Extract long sequences.

Write out lines of FASTQ and FASTA files that exceed an argument-specified
length.

% scripts/extract-long-sequences.py [-h] [-o OUTPUT] [-l LENGTH]
                                 input_filenames [input_filenames ...]

Use '-h' for parameter help.
"""
from __future__ import print_function
import argparse
import screed
import sys
from khmer.utils import write_record


def get_parser():
    parser = argparse.ArgumentParser(
        description='Extract FASTQ or FASTA sequences longer than'
        ' specified length (default: 200 bp).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_filenames', help='Input FAST[AQ]'
                        ' sequence filename.', nargs='+')
    parser.add_argument('-o', '--output', help='The name of the output'
                        ' sequence file.', default="/dev/stdout")
    parser.add_argument('-l', '--length', help='The minimum length of'
                        ' the sequence file.',
                        type=int, default=200)
    return parser


def main():
    args = get_parser().parse_args()
    outfp = open(args.output, 'w')
    for filename in args.input_filenames:
        for record in screed.open(filename, parse_description=False):
            if len(record['sequence']) >= args.length:
                write_record(record, outfp)
    print('wrote to: ' + args.output, file=sys.stderr)

if __name__ == '__main__':
    main()
