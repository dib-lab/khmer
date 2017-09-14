#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2012-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org

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
                        default='-')
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
