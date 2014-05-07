#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed
import os
import khmer

K = 32
N_HT = 4
PARTITION_SIZE_LIMIT = 5

###


def main():
    ht_filename = sys.argv[1]
    contig_filename = sys.argv[2]

    print>>sys.stderr, 'loading ht from', ht_filename
    ht = khmer.new_counting_hash(K, 1, N_HT)
    ht.load(ht_filename)

    partition_counts = {}

    for record in screed.open(contig_filename):
        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'G')

        a, b, c = ht.get_median_count(seq)

        partition = record.name.strip().split()[-1]

        x = partition_counts.get(partition, [])
        x.append(a)
        partition_counts[partition] = x

    for k, x in partition_counts.iteritems():
        if len(x) < PARTITION_SIZE_LIMIT:
            continue

        fp = open('partition%s.counts' % k, 'w')
        for i in x:
            fp.write("%s\n" % i)
        fp.close()

if __name__ == '__main__':
    main()
