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
Annotate sequences with partition numbers.

% python scripts/annotate-partitions.py <pmap_file> <file1> [ <file2> ... ]

Partition-annotated sequences will be in <fileN>.part.

Use '-h' for parameter help.
"""

import os
import textwrap
import sys
from khmer import __version__, Nodegraph
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import (sanitize_help, KhmerArgumentParser)

DEFAULT_K = 32


def get_parser():
    epilog = """\
    Load in a partitionmap (generally produced by :program:`partition-graph.py`
    or :program:`merge-partitions.py`) and annotate the sequences in the given
    files with their partition IDs. Use :program:`extract-partitions.py` to
    extract sequences into separate group files.

    Example (results will be in ``random-20-a.fa.part``)::

        load-graph.py -k 20 example tests/test-data/random-20-a.fa
        partition-graph.py example
        merge-partitions.py -k 20 example
        annotate-partitions.py -k 20 example tests/test-data/random-20-a.fa
    """
    parser = KhmerArgumentParser(
        description="Annotate sequences with partition IDs.",
        epilog=textwrap.dedent(epilog))

    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('graphbase', help='basename for input and output '
                        'files')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+', help='input FAST[AQ] sequences to '
                        'annotate.')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()

    ksize = args.ksize
    filenames = args.input_filenames
    nodegraph = Nodegraph(ksize, 1, 1)

    partitionmap_file = args.graphbase + '.pmap.merged'

    check_input_files(partitionmap_file, args.force)
    for _ in filenames:
        check_input_files(_, args.force)

    check_space(filenames, args.force)

    print('loading partition map from:', partitionmap_file, file=sys.stderr)
    nodegraph.load_partitionmap(partitionmap_file)

    for infile in filenames:
        print('outputting partitions for', infile, file=sys.stderr)
        outfile = os.path.basename(infile) + '.part'
        part_count = nodegraph.output_partitions(infile, outfile)
        print('output %d partitions for %s' % (
            part_count, infile), file=sys.stderr)
        print('partitions are in', outfile, file=sys.stderr)


if __name__ == '__main__':
    main()
