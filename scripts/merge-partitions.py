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
Merge multiple pmap files into a single one.

% python scripts/merge-partitions.py <base>

Load <base>.subset.*.pmap and merge into a single pmap file.  Final
merged pmap file will be in <base>.pmap.merged.
"""

import glob
import os
import textwrap
import khmer
import sys
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import sanitize_help, KhmerArgumentParser

DEFAULT_K = 32


def get_parser():
    epilog = """\
    Take the ``${graphbase}.subset.#.pmap`` files and merge them all into a
    single ``${graphbase}.pmap.merged`` file for
    :program:`annotate-partitions.py` to use.
    """
    parser = KhmerArgumentParser(
        description="Merge partition map '.pmap' files.",
        epilog=textwrap.dedent(epilog),
        citations=['graph'])
    parser.add_argument('-k', '--ksize', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('--keep-subsets', dest='remove_subsets',
                        default=True, action='store_false',
                        help='Keep individual subsets (default: False)')
    parser.add_argument('graphbase', help='basename for input and output '
                        'files')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()

    output_file = args.graphbase + '.pmap.merged'
    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print('loading %d pmap files (first one: %s)' %
          (len(pmap_files), pmap_files[0]), file=sys.stderr)

    ksize = args.ksize
    nodegraph = khmer.Nodegraph(ksize, 1, 1)

    for _ in pmap_files:
        check_input_files(_, args.force)

    check_space(pmap_files, args.force)

    for pmap_file in pmap_files:
        print('merging', pmap_file, file=sys.stderr)
        nodegraph.merge_subset_from_disk(pmap_file)

    print('saving merged to', output_file, file=sys.stderr)
    nodegraph.save_partitionmap(output_file)

    if args.remove_subsets:
        print('removing pmap files', file=sys.stderr)
        for pmap_file in pmap_files:
            os.unlink(pmap_file)


if __name__ == '__main__':
    main()
