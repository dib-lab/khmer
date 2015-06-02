#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import sys
import khmer
import screed
import os

K = 20
HASHTABLE_SIZE = int(4e9)
N_HT = 4

UNIQUE_LEN = 100
UNIQUE_F = 0.9


def main():
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    uniq2 = open(os.path.basename(sys.argv[2]) + '.uniq', 'w')

    kh = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
    for n, record in enumerate(screed.open(filename1)):
        if n % 10000 == 0:
            print('...', filename1, n)
        seq = record.sequence.upper().replace('N', 'A')
        kh.consume(seq)

    path_n = 0
    for n, record in enumerate(screed.open(filename2)):
        if n % 10000 == 0:
            print('...', filename2, n)
        seq = record.sequence.upper().replace('N', 'A')
        paths = kh.extract_unique_paths(seq, UNIQUE_LEN, UNIQUE_F)
        kh.consume(seq)

        for path in paths:
            path_n += 1
            print('>%s from:%s\n%s' % (path_n, record.name, path), file=uniq2)


if __name__ == '__main__':
    main()
