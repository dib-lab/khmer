#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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
Sequence trimming using stoptags.

Trim sequences at k-mers in the given stoptags file.  Output sequences
will be placed in 'infile.stopfilt'.

% python scripts/filter-stoptags.py <stoptags> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import os
import textwrap
import sys
from khmer import Nodegraph
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import sanitize_help, KhmerArgumentParser

# @CTB K should be loaded from file...
DEFAULT_K = 32


def get_parser():
    epilog = """\
    Load stoptags in from the given `.stoptags` file and use them to trim
    or remove the sequences in `<file1-N>`.  Trimmed sequences will be placed
    in `<fileN>.stopfilt`.
    """
    parser = KhmerArgumentParser(
        description="Trim sequences at stoptags.",
        epilog=textwrap.dedent(epilog), citations=['graph'])
    parser.add_argument('-k', '--ksize', default=DEFAULT_K, type=int,
                        help='k-mer size')
    parser.add_argument('stoptags_file', metavar='input_stoptags_filename')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()
    stoptags = args.stoptags_file
    infiles = args.input_filenames

    for _ in infiles:
        check_input_files(_, args.force)

    check_space(infiles, args.force)

    print('loading stop tags, with K', args.ksize, file=sys.stderr)
    nodegraph = Nodegraph(args.ksize, 1, 1)
    nodegraph.load_stop_tags(stoptags)

    def process_fn(record):
        name = record.name
        seq = record.sequence
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = nodegraph.trim_on_stoptags(seq)

        if trim_at >= args.ksize:
            return name, trim_seq

        return None, None

    # the filtering loop
    for infile in infiles:
        print('filtering', infile, file=sys.stderr)
        outfile = os.path.basename(infile) + '.stopfilt'

        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), outfp)

        print('output in', outfile, file=sys.stderr)


if __name__ == '__main__':
    main()
