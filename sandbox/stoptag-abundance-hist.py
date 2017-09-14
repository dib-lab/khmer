#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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
import os
try:
    import matplotlib
    matplotlib.use('Agg')
    from pylab import *
except ImportError:
    pass

def main():

    hashfile = sys.argv[1]
    filename = sys.argv[2]
    figure = sys.argv[3]

    ht = khmer.load_countgraph(hashfile)

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
