#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Use '-h' for parameter help.
"""
import os
import khmer
import sys
from khmer.khmer_args import build_counting_args
from khmer.utils import iter_kmers
import screed

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2

DEBUG = False

def write_seq(fp, name, seq, acc):
    if acc:
        fp.write('@{name}\n{seq}\n+\n{acc}\n'.format(
            name=name, seq=seq, acc=acc))
    else:
        fp.write('>{name}\n{seq}\n'.format(
            name=name, seq=seq))


def main():
    parser = build_counting_args()
    parser.add_argument('--cutoff', '-C', dest='cutoff',
                        default=DEFAULT_CUTOFF, type=int,
                        help="Trim at k-mers below this abundance.")

    parser.add_argument('-l', '--limit', type=int,
                        dest='trusted_limit',
                        help='position where we stop assuming k-mers are trusted,\
                                relative to the end of the read')
    parser.add_argument('--normalize-to', '-Z', type=int, dest='normalize_to',
                        help='base variable-coverage cutoff on this median'
                        ' k-mer abundance',
                        default=DEFAULT_NORMALIZE_LIMIT)
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()

    infiles = args.input_files

    C = args.cutoff
    Z = args.normalize_to

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    if args.trusted_limit is None:
        limit = K * 2
    else:
        limit = args.trusted_limit

    print >>sys.stderr, 'building hashtable...'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    print >>sys.stderr, "K:", K

    for infile in infiles:
        print >>sys.stderr, 'filtering', infile
        outfile = os.path.basename(infile) + '.kmerfilt'
        outbadkmers = os.path.basename(infile) + '.badkmers.fa'
        outgoodkmers = os.path.basename(infile) + '.goodkmers.fa'
        r_outfp = open(outfile, 'wb')
        k_bad_outfp = open(outbadkmers, 'wb')
        k_good_outfp = open(outgoodkmers, 'wb')
        for n, record in enumerate(screed.open(infile)):
            if n % 250000 == 0 and not args.quiet:
                print >>sys.stderr, 'processed {n} of {f}...'.format(n=n, f=infile)
            name = record['name']
            seq = record['sequence'].upper()

            if hasattr(record, 'accuracy'):
                acc = record['accuracy']
            else:
                acc = None
            if DEBUG:
                print 'checking', name, 'length', len(seq)

            if 'N' in seq:
                continue

            med, _, _ = ht.get_median_count(seq)
            if med < Z:
                if DEBUG:
                    print 'wrote ', name, 'to reads, med=', med
                write_seq(r_outfp, name, seq, acc)
                ht.consume(seq)
            else:
                if DEBUG:
                    print name, '>= Z, med=', med
                pos = max(0, len(seq) - limit + 1)
                if len(seq) < K:
                    if DEBUG:
                        print 'passed length filter'
                    continue 
                trimat = -1
                while (pos < (len(seq) - K + 1)):
                    if DEBUG:
                        print 'checking', name, 'pos', pos
                    kmer = seq[pos:pos+K]
                    #print kmer
                    kmer_c = ht.get(kmer)
                    if DEBUG:
                        print kmer, 'at pos', pos, 'C=', kmer_c
                    if kmer_c < C:
                        if DEBUG:
                            print trimat
                        if trimat < 0:
                            trimat = pos
                        write_seq(k_bad_outfp, '{name}-{pos}:{c}'.format(name=name, pos=pos, c=kmer_c), kmer, None)
                    else:
                        write_seq(k_good_outfp, '{name}-{pos}:{c}'.format(name=name, pos=pos, c=kmer_c), kmer, None)
                    pos += 1
                seq = seq[:trimat]
                if acc:
                    acc = acc[:trimat]
                write_seq(r_outfp, name, seq, acc)

        print >>sys.stderr, 'output in', outfile

if __name__ == '__main__':
    main()
