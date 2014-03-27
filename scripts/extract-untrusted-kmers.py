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

def print_dbg(*args):
    if DEBUG:
        print >>sys.stderr, ' '.join([str(a) for a in args])

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

    print >>sys.stderr, 'Using cutoff {c}, filter {z}, on last {l} bases'.format(
        c=C, z=Z, l=limit)

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

        n_bad = 0
        n_short = 0
        n_trimmed = 0
        total_out = 0

        for n, record in enumerate(screed.open(infile)):
            if n % 250000 == 0 and not args.quiet:
                print >>sys.stderr, 'processed {n} of {f}...'.format(n=n, f=infile)
            name = record['name']
            seq = record['sequence'].upper()

            if hasattr(record, 'accuracy'):
                acc = record['accuracy']
            else:
                acc = None
            print_dbg('checking', name, 'length', len(seq))

            if 'N' in seq or len(seq) < K:
                n_bad += 1
                continue

            med, _, _ = ht.get_median_count(seq)
            if med < Z:
                total_out += 1
                print_dbg('wrote ', name, 'to reads, med=', med)
                write_seq(r_outfp, name, seq, acc)
                ht.consume(seq)
            else:
                print_dbg(name, '>= Z, med=', med)
                pos = max(0, len(seq) - limit + 1)
                trimat = -1
                while (pos < (len(seq) - K + 1)):
                    print_dbg('checking', name, 'pos', pos)
                    kmer = seq[pos:pos+K]
                    kmer_c = ht.get(kmer)
                    print_dbg(kmer, 'at pos', pos, 'C=', kmer_c)
                    if kmer_c < C:
                        print_dbg(trimat)
                        if trimat < 0:
                            print_dbg(kmer, len(kmer), kmer_c)
                            #print kmer, len(kmer), kmer_c
                            trimat = pos
                        write_seq(k_bad_outfp, '{name}-{pos}:{c}'.format(name=name, pos=pos, c=kmer_c), kmer, None)
                    else:
                        write_seq(k_good_outfp, '{name}-{pos}:{c}'.format(name=name, pos=pos, c=kmer_c), kmer, None)
                    pos += 1
                if trimat >= K:
                    n_trimmed += 1
                    seq = seq[:trimat]
                    if acc:
                        acc = acc[:trimat]
                elif trimat > 0:
                    n_short += 1
                    continue

                total_out += 1
                write_seq(r_outfp, name, seq, acc)
        
        print total_out, 'passed with', n_trimmed, 'trimmed,', n_bad, 'failed filter', n_short, 'failed filter after trimming'
        print >>sys.stderr, 'output in', outfile
        

if __name__ == '__main__':
    main()
