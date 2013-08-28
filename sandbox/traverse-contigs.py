#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
"""

import khmer
import sys
import os.path

import sys
import screed
import khmer
from khmer.hashbits_args import build_construct_args, DEFAULT_MIN_HASHSIZE


def main():
    parser = build_construct_args()
    parser.add_argument('--build-tagset', '-t', default=True,
                        action='store_false',
                        help='Construct tagset while loading sequences')
    parser.add_argument('input_contigs')
    parser.add_argument('input_reads', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % \
            args.min_hashsize
        print >>sys.stderr, ''
        print >>sys.stderr, 'Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize / 8)' % (
            args.n_hashes * args.min_hashsize / 8.)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    base = args.input_contigs
    filenames = args.input_reads

    ###

    print 'making hashtable'
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    for filename in filenames:
        print 'consuming input', filename
        ht.consume_fasta(filename)

    ###

    print 'reading contigs from', args.input_contigs

    fp = open(os.path.basename(args.input_contigs) + '.unfound', 'w')

    for contig in screed.open(base):
        seq = contig.sequence

        n_not_found = 0
        found = True
        for start in range(0, len(seq) - K):
            kmer = seq[start:start + K]
            if not ht.get(kmer):
                n_not_found += 1
                found = False

        if not found:
            fp.write('>%s %d %d\n%s\n' % (contig.name,
                     n_not_found, len(contig.sequence), contig.sequence))

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
