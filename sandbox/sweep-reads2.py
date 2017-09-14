#!/usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2012-2015, Michigan State University.
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
"""
Use a set of query reads to sweep out overlapping reads from another file.

% python scripts/sweep-reads2.py <query reads> <search reads>

Results end up in <search reads>.sweep2.

Use '-h' for parameter help.
"""

import sys
import khmer
import os.path
import screed
from khmer import khmer_args
from khmer.khmer_args import (build_nodegraph_args, DEFAULT_MAX_TABLESIZE)
from khmer.utils import broken_paired_reader, write_record


def main():
    parser = build_nodegraph_args()
    parser.add_argument('-o', '--outfile',
                        help='output file; default is "infile".sweep2')
    parser.add_argument('-q', '--quiet')
    parser.add_argument('input_filename')
    parser.add_argument('read_filename')

    args = parser.parse_args()

    inp = args.input_filename
    readsfile = args.read_filename

    outfile = os.path.basename(readsfile) + '.sweep2'
    if args.outfile:
        outfile = args.outfile
    outfp = open(outfile, 'w')

    # create a nodegraph data structure
    ht = khmer_args.create_countgraph(args)

    # load contigs, connect into N partitions
    print('loading input reads from', inp)
    ht.consume_seqfile(inp)

    print('starting sweep.')

    m = 0
    K = ht.ksize()
    instream = screed.open(readsfile)
    for n, is_pair, read1, read2 in broken_paired_reader(instream):
        if n % 10000 == 0:
            print('...', n, m)

        if is_pair:
            count1 = ht.get_median_count(read1.sequence)[0]
            count2 = ht.get_median_count(read2.sequence)[0]
            if count1 or count2:
                m += 1
                write_record_pair(read1, read2, outfp)
        else:
            count = ht.get_median_count(read1.sequence)[0]
            if count:
                m += 1
                write_record(read1, outfp)

if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
