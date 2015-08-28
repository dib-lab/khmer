#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import screed
import sys
import textwrap
import argparse

from khmer import __version__
from khmer.kfile import (add_output_compression_type, get_file_writer,
                         is_block, describe_file_handle)
from khmer.utils import write_record

from khmer.khmer_args import (sanitize_help, ComboFormatter, info,
                              _VersionStdErrAction)

CUTOFF = 200

def get_parser():
    parser = argparse.ArgumentParser(
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).',
        formatter_class=ComboFormatter)

    parser.add_argument('tag', help='The name of the tag'
                        ' for example: assembly')
    parser.add_argument('input_filename', nargs='+',
                        help="The names of one or more groups "
                        "sequence files.")   
    parser.add_argument('-o', '--output', metavar="filename",
                        type=argparse.FileType('wb'),
                        help='The name of the output'
                        ' FASTA sequence file.',
                        default=sys.stdout)
    parser.add_argument('--version', action=_VersionStdErrAction,
                        version='khmer {v}'.format(v=__version__))
    add_output_compression_type(parser)
    return parser



def main():
    n = 0
    args = sanitize_help(get_parser()).parse_args()
    prefix = args.tag
    print (args.input_filename)
    for filename in args.input_filename:
       for record in screed.open(filename):
            if len(record.sequence) >= CUTOFF:
                n += 1
                print('>%s.%s %s' % (prefix, n, record.name))
                x = "\n".join(textwrap.wrap(record.sequence, 80))
                print (x)
   
if __name__ == '__main__':
    main()
