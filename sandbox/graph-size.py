#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
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
    ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    print('eating fa', infile)
    total_reads, n_consumed = ht.consume_fasta(infile)
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
