#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
'''
assemstats.py

Uses screed to calculate various assembly statistics.

You can obtain screed by running
   pip install screed
'''

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

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
