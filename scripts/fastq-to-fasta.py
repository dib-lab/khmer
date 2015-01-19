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
from oxli import fq2fa


def get_parser():
    parser = argparse.ArgumentParser(
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    fq2fa.add_args(parser)
    return parser


def main():
    args = get_parser().parse_args()
    fq2fa.do_fastq_to_fasta(args.input_sequence, args.output, args.n_keep)

if __name__ == '__main__':
    main()
