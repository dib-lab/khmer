#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
'''
assemstats.py

Uses screed to calculate various assembly statistics.

You can obtain screed at github by running
   git clone git://github.com/acr/screed.git

Then, install by running
   python setup.py install
in the newly created screed directory.

Once completed, you should be able to run this script as is.

Author: Jason Pell (pelljaso@cse.msu.edu)
'''
from __future__ import division
from __future__ import print_function

import screed
import sys
import glob
import os

def trimLens(lens, minLen):
    '''
    Eliminates any reads below a certain threshold.
    Function assumes that input list lens is sorted smallest to largest.
    '''

    index = 0

    for i in range(len(lens)):
        if lens[i] < minLen:
            index += 1
        else:
            break

    return lens[index:len(lens)]


def getLens(filename):
    '''
    Parses FASTA file using screed to create a sorted list of contig lengths.
    '''
    lens = []

    fa_instance = screed.open(filename)
    for record in fa_instance:
        lens.append(len(record['sequence']))

    return sorted(lens)


def calcNXX(lens, percent):
    '''
    Calculates any NXX (e.g. N50, N90) statistic.
    '''

    lenSum = sum(lens)
    threshold = (float(percent) / 100) * lenSum

    runningSum = 0
    nxx = 0
    nxxLen = 0

    for i in range(len(lens) - 1, -1, -1):
        myLen = lens[i]
        nxx += 1
        runningSum += myLen

        if runningSum >= threshold:
            nxxLen = myLen
            break

    return nxx, nxxLen


def main():
    '''
    Outputs assembly statistics for provided FASTA files.
    '''
    totalN = 0
    totalSum = 0

    if len(sys.argv) < 3:
        print("Usage: python assemstats.py <min contig length> [ FASTA files ]")
        return

    try:
        minLen = int(sys.argv[1])
    except ValueError:
        print("Minimum contig length must be an integer.")
        return

    print('** cutoff:', minLen)
    print("N\tsum\tmax\tfilename")

    for filename in sys.argv[2:]:
        if not os.path.exists(filename):
            print("WARNING: file %s does not exist." % filename, file=sys.stderr)
            continue

        lens = getLens(filename)
        trimmedLens = trimLens(lens, minLen)

        if trimmedLens:
            statTrimmedN = len(trimmedLens)
            statSum = sum(trimmedLens)
            statMax = max(trimmedLens)
        else:
            statTrimmedN = 0
            statSum = 0
            statMax = 0

        totalN += statTrimmedN
        totalSum += statSum

        print("%d\t%d\t%d\t%s" % (statTrimmedN, statSum, statMax, filename))

    if len(sys.argv) > 3 and totalN:
        print('--')
        print('TOTAL: %g in %d contigs (mean size %d)' % (
            totalSum, totalN, totalSum / totalN + .5))

main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
