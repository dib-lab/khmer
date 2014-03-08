#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

K = 32
HT_SIZE = 7e9
N_HT = 4

WORKER_THREADS = 8
GROUPSIZE = 100

###


def main():
    print '-- settings:'
    print 'K', K
    print 'N THREADS', WORKER_THREADS
    print '--'

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    for filename in sys.argv[1:]:
        print 'consuming input', filename
        ht.consume_fasta(filename)

    def process_fn(record, ht=ht):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        if len(seq) < K:
            return None, None

        if ht.get_min_count(seq) < 2:
            return None, None

        return name, seq

    for filename in sys.argv[1:]:
        print '***', filename
        outfile = os.path.basename(filename) + '.f2'
        if os.path.exists(outfile):
            print 'SKIPPING', outfile, ' -- already exists'
            continue

        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)
        tsp.start(verbose_fasta_iter(filename), outfp)

if __name__ == '__main__':
    main()
