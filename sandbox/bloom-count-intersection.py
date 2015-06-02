from __future__ import print_function
#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# using bloom filter to count intersection

import khmer
import sys
import screed
from screed.fasta import fasta_iter


def main():
    filename = sys.argv[1]
    K = int(sys.argv[2])  # size of kmer
    HT_SIZE = int(sys.argv[3])  # size of hashtable
    N_HT = int(sys.argv[4])  # number of hashtables

    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    n_unique = 0
    for n, record in enumerate(fasta_iter(open(filename))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            if (not ht.get(kmer)):
                n_unique += 1
            ht.count(kmer)
    print(filename, 'has been consumed.')
    print('# of unique kmers:', n_unique)
    print('# of occupied bin:', ht.n_occupied())

    filename2 = sys.argv[5]
    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT)
    n_unique = 0
    n_overlap = 0
    for n, record in enumerate(fasta_iter(open(filename2))):
        sequence = record['sequence']
        seq_len = len(sequence)
        for n in range(0, seq_len + 1 - K):
            kmer = sequence[n:n + K]
            if (not ht2.get(kmer)):
                n_unique += 1
                if (ht.get(kmer)):
                    n_overlap += 1
            ht2.count(kmer)

    print(filename2, 'has been consumed.')
    print('# of unique kmers:', n_unique)
    print('# of occupied bin:', ht2.n_occupied())

    print(n_overlap, 'unique kmers appears in both ', filename, ' and ', filename2)


if __name__ == '__main__':
    main()
