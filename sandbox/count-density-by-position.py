#! /usr/bin/env
import sys
import screed.fasta
import os
import khmer
import argparse

###

DEFAULT_MAX_DENSITY = 2000
DEFAULT_RADIUS = 10
DEFAULT_NUM_READS = int(1e6)

MAX_READ_LENGTH = 200

###

def main():
    parser = argparse.ArgumentParser(description=\
                 'Calculate local graph density by starting position in read')

    parser.add_argument('-M', '--max-density', type=int, dest='max_density',
                        default=DEFAULT_MAX_DENSITY)
    parser.add_argument('-r', '--radius', type=int, dest='radius',
                        default=DEFAULT_RADIUS)
    parser.add_argument('-n', '--num-reads', type=int, dest='num_reads',
                        default=DEFAULT_NUM_READS,
                        help="Number of reads to use to calculate average")

    parser.add_argument('htfile')
    parser.add_argument('input')
    parser.add_argument('output')

    args = parser.parse_args()

    htfile = args.htfile
    infile = args.input
    outfile = args.output

    print 'saving to:', outfile

    ht = khmer.new_hashbits(1, 1, 1)
    print 'loading hashtable', htfile
    ht.load(htfile)
    K = ht.ksize()
    RADIUS = args.radius

    print 'K is', K
    print 'RADIUS is', RADIUS

    hist = [0.0] * MAX_READ_LENGTH
    histcount = [0] * MAX_READ_LENGTH

    for n, record in enumerate(screed.fasta.fasta_iter(open(infile))):
        seq = record['sequence']

        for pos in range(0, len(seq) - K + 1):
            density = ht.count_kmers_within_radius(seq[pos:pos + K], RADIUS,
                                                   args.max_density)
            density /= float(RADIUS)
            hist[pos] += density
            histcount[pos] += 1

        if n % 1000 == 0:
            print '... saving', n

            outfp = open(outfile, 'w')
            for i in range(len(hist)):
                if histcount[i]:
                    print >>outfp, i, hist[i], histcount[i], \
                        hist[i] / histcount[i]

            if n >= args.num_reads:
                break

if __name__ == '__main__':
    main()
