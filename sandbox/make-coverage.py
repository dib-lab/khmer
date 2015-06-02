#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.txt.
# Contact: khmer-project@idyll.org
#

from __future__ import print_function

import screed

import sys

def main():
    dbfile = sys.argv[1]
    mapfile = sys.argv[2]

    lengths = {}
    for n, record in enumerate(screed.open(dbfile)):
        if n % 100000 == 0:
            print('...', n)
        lengths[record.name] = len(record.sequence)

    sums = {}
    for n, line in enumerate(open(mapfile)):
        if n % 100000 == 0:
            print('... 2x', n)
        x = line.split('\t')
        name = x[2]
        readlen = len(x[4])
        sums[name] = sums.get(name, 0) + 1

    mapped_reads = n

    rpkms = {}
    for k in sums:
        rpkms[k] = sums[k] * (1000. / float(lengths[k])) * \
            float(mapped_reads) / 1e6

    outfp = open(dbfile + '.cov', 'w')
    for n, record in enumerate(screed.open(dbfile)):
        if n % 100000 == 0:
            print('...', n)

        print(">%s[cov=%d]\n%s" % (record.name,
                                   rpkms.get(record.name, 0),
                                   record.sequence),
              file=outfp)

if __name__ == '__main__':
        main()
