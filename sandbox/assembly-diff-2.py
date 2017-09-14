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

    kh = khmer.Nodegraph(K, HASHTABLE_SIZE, N_HT)
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
