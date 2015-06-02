#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import sys
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_fasta_iter

WORKER_THREADS = 8
GROUPSIZE = 100

CUTOFF = 50

###


def main():
    counting_ht = sys.argv[1]
    infiles = sys.argv[2:]

    print('file with ht: %s' % counting_ht)
    print('-- settings:')
    print('N THREADS', WORKER_THREADS)
    print('--')

    print('making hashtable')
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    for infile in infiles:
        print('filtering', infile)
        outfile = os.path.basename(infile) + '.below'

        outfp = open(outfile, 'w')

        def process_fn(record, ht=ht):
            name = record['name']
            seq = record['sequence']
            if 'N' in seq:
                return None, None

            trim_seq, trim_at = ht.trim_below_abundance(seq, CUTOFF)

            if trim_at >= K:
                return name, trim_seq

            return None, None

        tsp = ThreadedSequenceProcessor(process_fn, WORKER_THREADS, GROUPSIZE)

        tsp.start(verbose_fasta_iter(infile), outfp)

if __name__ == '__main__':
    main()
