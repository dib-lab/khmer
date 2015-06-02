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
import os
try:
    from pylab import *
except ImportError:
    pass

def main():

    hashfile = sys.argv[1]
    filename = sys.argv[2]
    figure = sys.argv[3]

    ht = khmer.load_counting_hash(hashfile)

    outabund = open(os.path.basename(filename) + '.counts', 'w')

    counts = []
    d = {}
    for sequence in open(sys.argv[2]):
        sequence = sequence.strip()

        count = ht.get(sequence)
        counts.append(count)
        d[count] = d.get(count, 0) + 1

        if count > 1000:
            print(sequence, count, file=outabund)

    outfp = open(figure + '.countshist', 'w')
    sofar = 0
    sofar_cumu = 0
    for k in sorted(d.keys()):
        sofar += d[k]
        sofar_cumu += k * d[k]
        print(k, d[k], sofar, sofar_cumu, file=outfp)

    hist(counts, normed=True, cumulative=True, bins=100, range=(1, 1000))
    savefig(figure)


if __name__ == '__main__':
    main()
