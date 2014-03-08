#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

MAX_SODD = 3

WORKER_THREADS = 8
GROUPSIZE = 100

###


def main():
    htfile = sys.argv[1]
    outfiles = sys.argv[2:]

    print 'loading hashbits'
    ht = khmer.load_hashbits(htfile)

    def process_fn(record, ht=ht):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)

        if trim_at >= ht.ksize():
            return name, trim_seq

        return None, None

    for filename in outfiles:
        outpath = os.path.basename(filename) + '.sodd'
        outfp = open(outpath, 'w')

        tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)
        tsp.start(verbose_fasta_iter(filename), outfp)

if __name__ == '__main__':
    main()
