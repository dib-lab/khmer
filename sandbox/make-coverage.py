#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
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
