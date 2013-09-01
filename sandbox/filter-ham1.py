#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

K = 32
WORKER_THREADS = 8
GROUPSIZE = 100

###


def main():
    counting_ht = sys.argv[1]
    infiles = sys.argv[2:]

    print 'file with ht: %s' % counting_ht
    print '-- settings:'
    print 'N THREADS', WORKER_THREADS
    print '--'

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, 1, 1)
    ht.load(counting_ht)

    for infile in infiles:
        print 'filtering', infile
        outfile = os.path.basename(infile) + '.ham1filt'

        outfp = open(outfile, 'w')

        def process_fn(record, ht=ht):
            name = record['name']
            seq = record['sequence']
            if 'N' in seq:
                return None, None

            for pos in range(len(seq) - K):
                kmer = seq[pos:pos + K]
                if ht.max_hamming1_count(kmer) > 2000:
                    trim_at = pos + K - 1
                    seq = seq[:trim_at]
                    break

            if len(seq) >= K:
                return name, seq

            return None, None

        tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

        tsp.start(verbose_fasta_iter(infile), outfp)

if __name__ == '__main__':
    main()
