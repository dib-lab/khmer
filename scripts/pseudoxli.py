#!/usr/bin/env python
#
# This file is a part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under the
# three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

"""
Single entry point script for khmer
"""

import argparse
import sys
from oxli import test
from oxli import fq2fa


def get_parser():
    """
    returns the parser object for the oxli subcommand handler
    """
    parser = argparse.ArgumentParser(
        description='Single entry point script for khmer',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(help='foo')

    # sub parser for fastq-to-fasta
    parser_fq2fa = subparsers.add_parser('fastq-to-fasta', help="Converts \
            FASTQ format(.fq) files to FASTA format(.fa)")

    parser_fq2fa.add_argument('input_sequence', help='The name of the input'
                              ' FASTQ sequence file.')
    parser_fq2fa.add_argument('-o', '--output', metavar="filename",
                              help='The name of the output'
                              ' FASTA sequence file.',
                              type=argparse.FileType('w'),
                              default=sys.stdout)
    parser_fq2fa.add_argument('-n', '--n_keep', default=False,
                              action='store_true', help='Option to drop reads containing \'N\'s \
            in input_sequence file.')
    parser_fq2fa.set_defaults(func=fastq_to_fasta)

    return parser


def fastq_to_fasta(args):
    """
    fastq_to_fasta subcommand handler function
    """
    print args
    print "Hello from fastq!"
    test.do_test()
    fq2fa.do_fastq_to_fasta(args)


def main():
    """
    main function; does the parsing and kicks off the subcommand
    """
    args = get_parser().parse_args()
    print "Hello from main!"
    args.func(args)

if __name__ == '__main__':
    main()
