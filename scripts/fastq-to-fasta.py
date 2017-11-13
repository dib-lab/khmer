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
Convert FASTQ files to FASTA format.

% python scripts/fastq-to-fasta.py [ -n -o ] <fastq_name>

Use '-h' for parameter help.
"""
import sys
import screed
from khmer import __version__
from khmer.kfile import (add_output_compression_type, get_file_writer,
                         describe_file_handle)
from khmer.utils import write_record
from khmer.khmer_args import sanitize_help, KhmerArgumentParser
from khmer.khmer_args import FileType as khFileType


def get_parser():
    parser = KhmerArgumentParser(
        description='Converts FASTQ format (.fq) files to FASTA format (.fa).')

    parser.add_argument('input_sequence', help='The name of the input'
                        ' FASTQ sequence file.')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=khFileType('wb'),
                        help='The name of the output'
                        ' FASTA sequence file.',
                        default=sys.stdout)
    parser.add_argument('-n', '--n_keep', default=False, action='store_true',
                        help='Option to keep reads containing \'N\'s in '
                             'input_sequence file. Default is to drop reads')
    add_output_compression_type(parser)
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()

    print('fastq from ', args.input_sequence, file=sys.stderr)
    outfp = get_file_writer(args.output, args.gzip, args.bzip)
    n_count = 0
    for n, record in enumerate(screed.open(args.input_sequence)):
        if n % 10000 == 0:
            print('...', n, file=sys.stderr)

        sequence = record['sequence']

        if 'N' in sequence:
            if not args.n_keep:
                n_count += 1
                continue

        del record['quality']
        write_record(record, outfp)

    print('\n' + 'lines from ' + args.input_sequence, file=sys.stderr)

    if not args.n_keep:
        print(str(n_count) + ' lines dropped.', file=sys.stderr)

    else:
        print('No lines dropped from file.', file=sys.stderr)

    print('Wrote output to', describe_file_handle(args.output),
          file=sys.stderr)


if __name__ == '__main__':
    main()
