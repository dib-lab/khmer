#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Bin reads by abundance, with some heuristics thrown in to reduce memory
usage.

CTB.
"""
import os.path
import khmer
import argparse
import screed

K=20
N=4
X=1e8

DEFAULT_MIN=20
DEFAULT_MAX=100
PCNT=50
DEFAULT_BINSIZE=1000

def output_single(r):
    if hasattr(r, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (r.name, r.sequence, r.accuracy)
    else:
        return ">%s\n%s\n" % (r.name, r.sequence)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("readfiles", nargs="+")
    parser.add_argument("-m", "--mincov", dest="mincov",
                        type=int, default=DEFAULT_MIN)
    parser.add_argument("-M", "--maxcov", dest="maxcov",
                        type=int, default=DEFAULT_MAX)
    parser.add_argument("-b", "--binsize", dest="binsize",
                        type=int, default=DEFAULT_BINSIZE)

    args = parser.parse_args()

    kh = khmer.new_counting_hash(K, X, N)
    elim = khmer.new_hashbits(K, X, N)

    for filename in args.readfiles:
        outp_index = 1
        outp_count = 0
        outputfilename = os.path.basename(filename) + '.bink.%d' % outp_index
        fp = open(outputfilename, 'w')

        n = 0
        n_low = 0
        n_med = 0
        n_high = 0
        n_elim = 0
        
        for n, record in enumerate(screed.open(filename)):
            if n % 25000 == 0:
                print '...', n
                
            seq = record.sequence.upper()
            if 'N' in seq:
                seq = seq.replace('N', 'G')

            if K <= len(seq):
                do_elim, _, _ = elim.get_median_count(seq)
                if do_elim:
                    n_elim += 1
                    continue
                
                a, _, _ = kh.get_median_count(seq)

                if a < args.mincov:              # low coverage? keep.
                    n_low += 1
                    kh.consume(seq)
                    do_output = True
                elif a >= args.maxcov:           # high coverage? discard.
                    n_high += 1
                    elim.consume(seq)
                    do_output = False
                else:                   # medium coverage? keep, sort of :)
                    n_med += 1
                    kh.consume_high_abund_kmers(seq, args.mincov)
                    do_output = True

            if do_output:
                outp_count += 1

         
                fp.write(output_single(record))

                if outp_count == args.binsize:
                    outp_index += 1
                    outputfilename = os.path.basename(filename) + \
                                     '.bink.%d' % outp_index
                    fp = open(outputfilename, 'w')
                    outp_count = 0

        print n + 1, "total reads"
        print n_low, "low abundance reads kept"
        print n_med, "med abundance reads kept"
        print n_elim, "high abundance reads eliminated"

if __name__ == '__main__':
    main()
