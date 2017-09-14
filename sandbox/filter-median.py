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
"""
Accept or discard sequences XXX, based on the given counting
hash table.  Output sequences will be placed in 'infile.medfilt'.

% python sandbox/filter-median.py <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

from khmer.khmer_args import build_counting_args

import random

###

DEFAULT_COVERAGE = 20


def main():
    parser = build_counting_args()
    parser.add_argument('--coverage', '-C', dest='coverage',
                        default=DEFAULT_COVERAGE, type=int)
    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print('file with ht: %s' % counting_ht)

    print('loading hashtable')
    ht = khmer.load_countgraph(counting_ht)
    K = ht.ksize()

    print("K:", K)

    # the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']

        med, avg, dev = ht.get_median_count(seq)

        if random.randint(1, med) > args.coverage:
            return None, None

        return name, seq

    # the filtering loop
    for infile in infiles:
        print('filtering', infile)
        outfile = os.path.basename(infile) + '.medfilt'
        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), outfp)

        print('output in', outfile)

if __name__ == '__main__':
    main()
