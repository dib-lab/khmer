#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
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
# pylint: disable=invalid-name,missing-docstring

"""
Extract long sequences.

Write out lines of FASTQ and FASTA files that exceed an argument-specified
length.

% scripts/extract-long-sequences.py [-h] [-o OUTPUT] [-l LENGTH]
                                 input_filenames [input_filenames ...]

Use '-h' for parameter help.
"""
import argparse
import screed
import textwrap
import sys
from khmer import __version__
from khmer.utils import write_record
from khmer.kfile import add_output_compression_type, get_file_writer
from khmer.khmer_args import sanitize_help, KhmerArgumentParser


def get_parser():
    epilog = """\
    Example::

        extract-long-sequences.py --length 10 tests/test-data/paired-mixed.fa
    """
    parser = KhmerArgumentParser(
        description='Extract FASTQ or FASTA sequences longer than'
        ' specified length (default: 200 bp).',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('input_filenames', help='Input FAST[AQ]'
                        ' sequence filename.', nargs='+')
    parser.add_argument('-o', '--output', help='The name of the output'
                        ' sequence file.', default=sys.stdout,
                        metavar='output', type=argparse.FileType('wb'))
    parser.add_argument('-l', '--length', help='The minimum length of'
                        ' the sequence file.',
                        type=int, default=200)
    add_output_compression_type(parser)
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()
    outfp = get_file_writer(args.output, args.gzip, args.bzip)
    for filename in args.input_filenames:
        for record in screed.open(filename):
            if len(record['sequence']) >= args.length:
                write_record(record, outfp)
    print('wrote to: ' + args.output.name, file=sys.stderr)


if __name__ == '__main__':
    main()
