#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
import khmer
import screed
import os
from screed.fasta import fasta_iter

K = 32
HASHTABLE_SIZE = 5e7
N_HT = 4

LENGTH_THRESHOLD = int(sys.argv[3])


def iterseq(x):
    for record in x:
        yield record.sequence


def slidingWindow(sequence, K_size):
    total_windows = (len(sequence) - K_size) + 1

    for i in range(0, total_windows, 1):
        kmer = sequence[i:i + K_size]
        if 'N' in kmer:
            continue
        yield sequence[i:i + K_size]


def count(contigs1, contigs2):
    ht = khmer.new_counting_hash(K, HASHTABLE_SIZE, N_HT)
    count = 0
    count2 = 0
    count3 = 0
    for n in contigs1:
        if len(n) >= LENGTH_THRESHOLD:
            kmer1 = slidingWindow(n, K)
            for x in kmer1:
                count += 1
                if ht.get(x):
                    continue
                ht.consume(x)

    for n in contigs2:
        if len(n) >= LENGTH_THRESHOLD:
            kmer2 = slidingWindow(n, K)
            for x in kmer2:
                count2 += 1
                if ht.get(x) > 0:
                    count3 += 1

    # 'count' is the total number of kmers in the first file
    # 'count2' is the total number of kmers in the second file
    # 'count3' is the total number of kmers shared between the two files.
    print count, count2, count3, "%.1f%%" % (count3 / float(count) * 100.)


def main(contig1, contig2):
    ht = count(iterseq(screed.open(contig1)),
               iterseq(screed.open(contig2)))


def output(ht, contigs, fp):
    for n, sequence in enumerate(contigs):
        count_list = []
        contig_kmers = slidingWindow(sequence, K)
        for contig_kmer in contig_kmers:
            count_kmer = int(ht.get(contig_kmer)) - 1
            count_list.append(count_kmer)

    for item in count_list:
        print >>fp, '%s' % item

    print 'Hashtable occupancy =', ht.n_occupied() / float(HASHTABLE_SIZE)
    return count_list

if __name__ == '__main__':

    contig_file1 = sys.argv[1]
    contig_file2 = sys.argv[2]

    main(contig_file1, contig_file2)
