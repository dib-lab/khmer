#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import division
from __future__ import print_function
import sys
import khmer
import screed
import os

K = 20
HASHTABLE_SIZE = int(2.5e8)
N_HT = 4

THRESHOLD = 0.9


def main():
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    uniq1 = open(os.path.basename(sys.argv[1]) + '.uniq', 'w')
    uniq2 = open(os.path.basename(sys.argv[2]) + '.uniq', 'w')
    paths = sys.argv[3]

    kh1 = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    kh1.consume_fasta(filename1)
    kh2 = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    kh2.consume_fasta(filename2)

    for record in screed.open(paths):
        n = 0
        n_present = 0

        path = record.sequence
        n = len(path) - K + 1
        for i in range(n):
            if kh1.get(path[i:i + K]):
                n_present += 1

        if n_present / float(n) >= THRESHOLD:
            present1 = True
        else:
            present1 = False

        n = 0
        n_present = 0

        path = record.sequence
        n = len(path) - K + 1
        for i in range(n):
            if kh2.get(path[i:i + K]):
                n_present += 1

        if n_present / float(n) >= THRESHOLD:
            present2 = True
        else:
            present2 = False

        if present1 and not present2:
            print('>%s\n%s' % (record.name, record.sequence), file=uniq1)
        elif present2 and not present1:
            print('>%s\n%s' % (record.name, record.sequence), file=uniq2)


if __name__ == '__main__':
    main()
