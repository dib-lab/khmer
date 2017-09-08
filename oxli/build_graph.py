#! /usr/bin/env python
#
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
# pylint: disable=missing-docstring
"""
Build a graph from the given sequences, save in <ptname>.

% python scripts/load-graph.py <ptname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""


import sys

import khmer
from khmer import khmer_args
from khmer.khmer_args import (report_on_config, info, add_threading_args,
                              calculate_graphsize)
from khmer.kfile import check_input_files
from khmer.kfile import check_space_for_graph
from oxli import functions as oxfuncs


def build_parser(parser):
    add_threading_args(parser)
    parser.add_argument('--no-build-tagset', '-n', default=False,
                        action='store_true', dest='no_build_tagset',
                        help='Do NOT construct tagset while loading sequences')
    parser.add_argument('output_filename',
                        metavar='output_nodegraph_filename', help='output'
                        ' k-mer nodegraph filename.')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+', help='input FAST[AQ] sequence filename')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main(args):
    graph_type = 'nodegraph'
    report_on_config(args, graphtype=graph_type)
    base = args.output_filename
    filenames = args.input_filenames

    for fname in args.input_filenames:
        check_input_files(fname, args.force)

    graphsize = calculate_graphsize(args, graph_type)
    space_needed = (args.n_tables * graphsize /
                    khmer._buckets_per_byte[graph_type])
    check_space_for_graph(args.output_filename, space_needed, args.force)

    print('Saving k-mer nodegraph to %s' % base, file=sys.stderr)
    print('Loading kmers from sequences in %s' %
          repr(filenames), file=sys.stderr)
    if args.no_build_tagset:
        print('We WILL NOT build the tagset.', file=sys.stderr)
    else:
        print('We WILL build the tagset (for partitioning/traversal).',
              file=sys.stderr)

    print('making nodegraph', file=sys.stderr)
    nodegraph = khmer_args.create_nodegraph(args)

    oxfuncs.build_graph(filenames, nodegraph, args.threads,
                        not args.no_build_tagset)

    print('Total number of unique k-mers: {0}'.format(
        nodegraph.n_unique_kmers()), file=sys.stderr)

    print('saving k-mer nodegraph in', base, file=sys.stderr)
    nodegraph.save(base)

    if not args.no_build_tagset:
        print('saving tagset in', base + '.tagset', file=sys.stderr)
        nodegraph.save_tagset(base + '.tagset')

    info_fp = open(base + '.info', 'w')
    info_fp.write('%d unique k-mers' % nodegraph.n_unique_kmers())

    fp_rate = \
        khmer.calc_expected_collisions(
            nodegraph, args.force, max_false_pos=.15)
    # 0.18 is ACTUAL MAX. Do not change.

    print('false positive rate estimated to be %1.3f' % fp_rate,
          file=sys.stderr)
    print('\nfalse positive rate estimated to be %1.3f' % fp_rate,
          file=info_fp)

    print('wrote to ' + base + '.info and ' + base, file=sys.stderr)
    if not args.no_build_tagset:
        print('and ' + base + '.tagset', file=sys.stderr)

    sys.exit(0)


if __name__ == '__main__':
    main(None)

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
