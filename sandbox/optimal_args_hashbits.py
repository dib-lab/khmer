#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2015, Michigan State University.
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
Estimate optimal arguments using nodegraph counting.

% python sandbox/optimal_args_nodegraph.py  <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import math
import threading

import khmer
from khmer.khmer_args import (report_on_config, info, add_threading_args,
                              build_nodegraph_args)
from khmer.kfile import check_input_files, check_space
from khmer.kfile import check_space
from khmer.khmer_args import graphsize_args_report


def get_parser():
    parser = build_nodegraph_args(descr="Load sequences into the compressible "
                                 "graph format plus optional tagset.")
    add_threading_args(parser)
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+', help='input FAST[AQ] sequence filename')
    return parser


def main():
    info('optimal_args_nodegraph.py', ['graph', 'SeqAn'])
    args = get_parser().parse_args()
    report_on_config(args, graphtype='nodegraph')


    filenames = args.input_filenames
    base = filenames[0]
    for _ in args.input_filenames:
        check_input_files(_, False)

    check_space(args.input_filenames, False)

    print('Counting kmers from sequences in %s' % repr(filenames),
          file=sys.stderr)

    htable = khmer.new_nodegraph(args.ksize, args.max_tablesize, args.n_tables)
    target_method = htable.consume_seqfile

    for _, filename in enumerate(filenames):
        rparser = khmer.ReadParser(filename)
        threads = []
        print('consuming input', filename, file=sys.stderr)
        for num in xrange(args.threads):
            cur_thread = threading.Thread(
                target=target_method, args=(rparser,))
            threads.append(cur_thread)
            cur_thread.start()

        for thread in threads:
            thread.join()
    unique_kmers = htable.n_unique_kmers()
    print('Total number of unique k-mers: {0}'.format(unique_kmers),
          file=sys.stderr)

    info_optimal = open(base + '.optimal_args', 'w')

    fp_rate = khmer.calc_expected_collisions(htable)
    print('fp rate estimated to be %1.3f' % fp_rate, file=sys.stderr)

    if fp_rate > 0.15:          # 0.18 is ACTUAL MAX. Do not change.
        print("**", file=sys.stderr)
        print("** ERROR: the graph structure is too small for this data set."
              "Increase table size/# tables.", file=sys.stderr)
        print("**", file=sys.stderr)
        if not False:
            sys.exit(1)

    to_print = graphsize_args_report(unique_kmers, fp_rate)
    
    print(to_print, file=info_optimal)
    
    print('optimal arguments were written to', base + '.optimal_args',
          file=sys.stderr)
    
if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
