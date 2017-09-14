#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2012-2015, Michigan State University.
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
import sys
import khmer
import argparse
import os
import screed


def main():
    parser = argparse.ArgumentParser(
        description="Output k-mer abundance distribution.")

    parser.add_argument('hashname')
    parser.add_argument('seqfile')
    parser.add_argument('histout')

    args = parser.parse_args()
    hashfile = args.hashname
    seqfile = args.seqfile
    histout = args.histout

    outfp = open(histout, 'w')

    print('hashtable from', hashfile)
    ht = khmer.load_countgraph(hashfile)

    hist = {}

    for i in range(65536):
        hist[i] = 0

    for n, record in enumerate(screed.open(seqfile)):
        if n > 0 and n % 100000 == 0:
            print('...', n)

        seq = record.sequence.replace('N', 'A')

        try:
            med, _, _ = ht.get_median_count(seq)
        except ValueError:
            continue

        hist[med] = hist[med] + 1

    histlist = list(hist.items())
    histlist.sort()

    maxk = max(hist.keys())
    sumk = sum(hist.values())

    sofar = 0
    for n, m in histlist:
        sofar += m
        percent = float(sofar) / sumk
        outfp.write('%d %d %d %.3f\n' % (n, m, sofar, percent))
    outfp.close()

if __name__ == '__main__':
    main()
