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
import khmer
import sys
import screed
import os.path
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

K = 32
HASHTABLE_SIZE = int(4e6)
THRESHOLD = 500
N_HT = 4
WORKER_THREADS = 5

###

GROUPSIZE = 100

###


def main():
    infile = sys.argv[1]
    outfile = os.path.basename(infile) + '.graphsize'
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

    print('input file to graphsize filter: %s' % infile)
    print('filtering to output:', outfile)
    print('-- settings:')
    print('K', K)
    print('HASHTABLE SIZE %g' % HASHTABLE_SIZE)
    print('N HASHTABLES %d' % N_HT)
    print('THRESHOLD', THRESHOLD)
    print('N THREADS', WORKER_THREADS)
    print('--')

    print('creating ht')
    ht = khmer.Nodegraph(K, HASHTABLE_SIZE, N_HT)
    print('eating fa', infile)
    total_reads, n_consumed = ht.consume_seqfile(infile)
    outfp = open(outfile, 'w')

    ###

    def process_fn(record, ht=ht):
        kmer = record['sequence'][:K]
        size = ht.calc_connected_graph_size(kmer, THRESHOLD)
        if size >= THRESHOLD:
            return record['name'], record['sequence']

        return None, None

    tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

    ###

    tsp.start(verbose_fasta_iter(infile), outfp)

if __name__ == '__main__':
    main()
