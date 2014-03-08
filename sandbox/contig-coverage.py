#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import khmer
import screed
import os
from screed.fasta import fasta_iter

K = 32
HASHTABLE_SIZE = 1e6
N_HT = 4

DATA_DIR = '/u/adina/test-reads/'

thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

datadir = os.path.join(thisdir, DATA_DIR)
datadir = os.path.abspath(datadir)


def slidingWindow(sequence, K_size):
    total_windows = (len(sequence) - K_size) + 1

    for i in range(0, total_windows, 1):
        yield sequence[i:i + K_size]


def main(contig_filename, read_filenames_list, output_filename):

    ht = khmer.new_counting_hash(K, HASHTABLE_SIZE, N_HT)

    '''consumes contig into hashtable'''
    for n, record in enumerate(fasta_iter(open(contig_filename))):
        sequence = record['sequence']
        contig_kmers = slidingWindow(sequence, K)
        for x in contig_kmers:
            if x.find('N') > 0:
                continue
            else:
                ht.consume(x)

    '''counts reads into hashtable abundance'''
    for each_file in read_filenames_list:
        read_file = open(each_file, 'r')
        for n1, record1 in enumerate(fasta_iter(read_file)):
            sequence = record1['sequence']
            read_kmers = slidingWindow(sequence, K)
            for kmer in read_kmers:
                if ht.get(kmer) > 0:
                    ht.count(kmer)
        read_file.close()

    '''retrieve abundances'''
    for n2, record2 in enumerate(fasta_iter(open(contig_filename))):
        contig_seq = record2['sequence']
        count_list = []
        contig_kmers = slidingWindow(contig_seq, K)
        for contig_kmer in contig_kmers:
            count_kmer = int(ht.get(contig_kmer)) - 1
            count_list.append(count_kmer)

    fp = open(output_filename, 'w')
    for item in count_list:
        print >>fp, '%s' % item

    print 'Hashtable occupancy =', ht.n_occupied() / float(HASHTABLE_SIZE)


def test():
    test_contig_file = os.path.join(datadir, 'sim-genome.fa')
    test_read_file = [os.path.join(datadir, 'fifty-reads-step3-10x.fa')]
    test_output_file = 'test.out'
    main(test_contig_file, test_read_file, test_output_file)


if '__name__==__main__':

    output_file = sys.argv[1]
    contig_file = sys.argv[2]
    reads_file = sys.argv[3:]

    # test()
    main(contig_file, reads_file, output_file)
