#! /usr/bin/env python
"""
Streaming error trimming based on digital normalization.

% python sandbox/trim-low-abund.py [ <data1> [ <data2> [ ... ] ] ]

Use -h for parameter help.
"""
import sys
import screed
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
import argparse

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2

DEFAULT_K = 32
DEFAULT_N_HT = 4
DEFAULT_MIN_HASHSIZE = 1e6


def main():
    parser = argparse.ArgumentParser(description='XXX')

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
    env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', DEFAULT_MIN_HASHSIZE)

    parser.add_argument('--ksize', '-k', type=int, dest='ksize',
                        default=env_ksize,
                        help='k-mer size to use')
    parser.add_argument('--n_hashes', '-N', type=int, dest='n_hashes',
                        default=env_n_hashes,
                        help='number of hash tables to use')
    parser.add_argument('--hashsize', '-x', type=float, dest='min_hashsize',
                        default=env_hashsize,
                        help='lower bound on hashsize to use')

    parser.add_argument('--cutoff', '-C', type=int, dest='abund_cutoff',
                        help='remove k-mers below this abundance',
                        default=DEFAULT_CUTOFF)

    parser.add_argument('--normalize-to', '-Z', type=int, dest='normalize_to',
                        help='base cutoff on median k-mer abundance of this',
                        default=DEFAULT_NORMALIZE_LIMIT)

    parser.add_argument('--mrna', '-m', dest='is_mrna',
                        help='treat as mRNAseq data',
                        default=True, action='store_true')

    parser.add_argument('--genome', '-g', dest='is_genome',
                        help='treat as genomic data (uniform coverage)',
                        default=False, action='store_true')

    parser.add_argument('--metagenomic', '-M',
                        dest='is_metagenomic',
                        help='treat as metagenomic data',
                        default=True, action='store_true')

    parser.add_argument('input_filenames', nargs='+')
    args = parser.parse_args()

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    CUTOFF = args.abund_cutoff
    NORMALIZE_LIMIT = args.normalize_to

    is_variable_abundance = True        # conservative
    if args.is_genome:
        is_variable_abundance = False

    errors = [0] * 1000

    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    save_pass2 = 0

    pass2list = []
    for filename in args.input_filenames:
        pass2filename = os.path.basename(filename) + '.pass2'
        trimfilename = os.path.basename(filename) + '.abundtrim'

        pass2list.append((pass2filename, trimfilename))

        pass2fp = open(pass2filename, 'w')
        trimfp = open(trimfilename, 'w')

        for n, read in enumerate(screed.open(filename)):
            if n % 10000 == 0:
                print '...', n, filename, save_pass2
            seq = read.sequence.replace('N', 'A')
            med, _, _ = ht.get_median_count(seq)

            if med < NORMALIZE_LIMIT:
                ht.consume(seq)
                pass2fp.write('>%s\n%s\n' % (read.name, read.sequence))
                save_pass2 += 1
            else:
                trim_seq, trim_at = ht.trim_on_abundance(seq, CUTOFF)
                if trim_at < len(seq):
                    errors[trim_at] += 1
                if trim_at >= K:
                    trimfp.write('>%s\n%s\n' % (read.name, trim_seq))

        pass2fp.close()
        trimfp.close()

        print 'saved %d of %d to pass2fp' % (save_pass2, n,)

    for pass2filename, trimfilename in pass2list:
        for n, read in enumerate(screed.open(pass2filename)):
            if n % 10000 == 0:
                print '... x 2', n, filename

            trimfp = open(trimfilename, 'a')

            seq = read.sequence.replace('N', 'A')
            med, _, _ = ht.get_median_count(seq)

            if med >= NORMALIZE_LIMIT or not is_variable_abundance:
                trim_seq, trim_at = ht.trim_on_abundance(seq, CUTOFF)
                if trim_at < len(seq):
                    errors[trim_at] += 1
                if trim_at >= K:
                    trimfp.write('>%s\n%s\n' % (read.name, trim_seq))
            else:
                trimfp.write('>%s\n%s\n' % (read.name, read.sequence))

    os.unlink(pass2filename)

    fp = open('err-profile.out', 'w')
    for pos, count in enumerate(errors):
        print >>fp, pos, count

if __name__ == '__main__':
    main()
