#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
# 
"""
Calculate the mismatch error profile for shotgun data, using a subset of
reads.

% scripts/calc-error-profile.py [ -o outfile ] <infile>

Reads FASTQ and FASTA input.  Output histogram is placed in
<infilename>.errhist OR in the filename passed to -o.
"""

import sys
import argparse
import khmer
import screed
import os.path

N_HT = 4
HASHSIZE = 1e8
K=20
C=5

MAX_SEQ_LEN = 65535
MAX_READS = 1e8

def exit_condition(n2_consumed, n_consumed):
    return (n2_consumed >= n_consumed or \
            n2_consumed > 1e5)

def main():
    parser = argparse.ArgumentParser(\
     "Calculate read error profile based on k-mer abundances of shotgun data.")

    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', dest='output_file',
                        type=argparse.FileType('w'), default=None)
    
    args = parser.parse_args()

    #
    # Figure out what the output filename is going to be
    #

    output_file = args.output_file
    if output_file:
        output_filename = output_file.name
    else:
        filename = args.filenames[0]
        output_filename = os.path.basename(filename) + '.errhist'
        output_file = open(output_filename, 'w')

    print 'Error histogram will be in %s' % output_filename

    ## Start!

    # build a small counting hash w/default parameters. In general there
    # should be no need to change these parameters.
    ht = khmer.new_counting_hash(K, HASHSIZE, N_HT)

    # initialize list to contain counts of errors by position
    positions = [0]*MAX_SEQ_LEN
    lengths = []                  # keep track of sequence lengths

    n_consumed = n2_consumed = 0
    bp_consumed = bp2_consumed = 0
    total = 0

    # run through all the files; pick out reads; once they saturate,
    # look for errors.
    total = 0
    for filename in args.filenames:
        print 'opening', filename

        for n, record in enumerate(screed.open(filename)):
            total += 1

            if total % 25000 == 0:
                print '...', total, n_consumed, n2_consumed

                # two exit conditions: first, have we hit our max reads limit?
                if total >= MAX_READS:
                    break

                # OR, alternatively, have we collected enough reads?
                if exit_condition(n2_consumed, n_consumed):
                    break

            # for each sequence, calculate its coverage:
            seq = record.sequence.replace('N', 'A')
            med, _, _ = ht.get_median_count(seq)

            # if the coverage is unsaturated, consume.
            if med < C:
                ht.consume(seq)
                n_consumed += 1
                bp_consumed += len(seq)
            else:
                # also consume & track up to 2C -- CTB consume only high abnd?
                if med < 2*C:
                    ht.consume(seq)
                    n2_consumed += 1
                    bp2_consumed += len(seq)
                    
                # for saturated data, find low-abund k-mers
                posns = ht.find_low_abund_kmers(seq, 2)
                lengths.append(len(seq))

                # track the positions => errors
                for p in posns:
                    positions[p] += 1

    # normalize for length
    lengths.sort()
    max_length = lengths[-1]

    length_count = [0]*max_length
    for j in range(max_length):
        length_count[j] = sum([ 1 for i in lengths if i >= j ])

    # find the last point for which we have any data
    last_zero = len(positions) - 1
    while positions[last_zero] == 0:
        last_zero -= 1

    # write!
    for n, i in enumerate(positions[:last_zero + 1]):
        print >>output_file, n, i, float(i) / float(length_count[n])

    output_file.close()

    print ''
    print 'total sequences:', total
    print 'n consumed:', n_consumed, n2_consumed
    print 'bp consumed:', bp_consumed, bp_consumed / float(C)

    if not exit_condition(n2_consumed, n_consumed):
        print >>sys.stderr, ""
        print >>sys.stderr, "** WARNING: not enough reads to get a good result"
        print >>sys.stderr, "** Is this high diversity sample / small subset?"
        sys.exit(-1)

if __name__ == '__main__':
    main()
