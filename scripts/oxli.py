#!/usr/bin/env python
#
# This file is a part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under the
# three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring

"""
Single entry point script for khmer
"""

import argparse

def get_parser():
    parser = argparse.ArgumentParser(
            description='Single entry point script for khmer',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_subparsers(help='foo')

    # sub parser for fastq-to-fasta
    parser_fq2fa = subparsers.add_parser('fastq-to-fasta', help='Converts FASTQ
    format (.fq) files to FASTA format (.fa)')

print("Hello, world; I'm oxli!")
