#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
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
    stoptags = sys.argv[1]
    infile = sys.argv[2]

    outfile = os.path.basename(infile) + '.stopkeep'
    if len(sys.argv) >= 4:
        outfile = sys.argv[3]

    print 'file with stop tags: %s' % stoptags
    print 'input file to filter: %s' % infile
    print 'filtering to output:', outfile
    print '-- settings:'
    print 'K', K
    print 'N THREADS', WORKER_THREADS
    print '--'

    print 'making hashtable'
    ht = khmer.new_hashbits(K, 1, 1)

    ht.load_stop_tags(stoptags)

    outfp = open(outfile, 'w')

    def process_fn(record, ht=ht):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = ht.trim_on_stoptags(seq)

        if trim_at < K:
            return name, seq

        seq = seq[trim_at:]
        if seq:
            return name, seq

        return None, None

    tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

    ###

    tsp.start(verbose_fasta_iter(infile), outfp)

if __name__ == '__main__':
    main()
