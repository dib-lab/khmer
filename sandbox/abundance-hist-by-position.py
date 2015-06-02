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

def main():
    n = 0
    countSum = [0] * 255
    countN = [0] * 255

    freqfile = sys.argv[1]

    print('opening .freq file:', freqfile, file=sys.stderr)
    fd = open(freqfile)
    for n, line in enumerate(fd):
        if n % 100000 == 0:
            print('...', n, file=sys.stderr)

        tok = line.split()

        for i in range(len(tok)):
            countSum[i] += int(tok[i])
            countN[i] += 1

    print('summarizing.', file=sys.stderr)
    y = [0.0] * len(countSum)

    for i in range(len(countSum)):
        if countN[i]:
            y[i] = float(countSum[i]) / float(countN[i])

    for n, i in enumerate(y):
        print(n, i)

if __name__ == '__main__':
    main()
