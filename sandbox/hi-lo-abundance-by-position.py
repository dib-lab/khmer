#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import sys
import os
import khmer


def write_dist(dist, fp):
    for n, i in enumerate(dist):
        fp.write('%d %d\n' % (n, i))


def main():
    hashfile = sys.argv[1]
    filename = sys.argv[2]
    outfile = os.path.basename(filename)

    print('loading kh file', hashfile)
    ht = khmer.load_counting_hash(hashfile)

    x = ht.fasta_count_kmers_by_position(filename, 100, 1)
    write_dist(x, open(outfile + '.pos.abund=1', 'w'))
    print('wrote', outfile + '.pos.abund=1')

    y = ht.fasta_count_kmers_by_position(filename, 100, 255)
    write_dist(y, open(outfile + '.pos.abund=255', 'w'))
    print('wrote', outfile + '.pos.abund=255')


if __name__ == '__main__':
    main()
