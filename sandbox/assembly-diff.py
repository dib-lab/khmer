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
HASHTABLE_SIZE = int(2.5e8)
N_HT = 4

THRESHOLD = 0.9


def main():
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    uniq1 = open(os.path.basename(sys.argv[1]) + '.uniq', 'w')
    uniq2 = open(os.path.basename(sys.argv[2]) + '.uniq', 'w')
    paths = sys.argv[3]

    kh1 = khmer.Nodegraph(K, HASHTABLE_SIZE, N_HT)
    kh1.consume_seqfile(filename1)
    kh2 = khmer.Nodegraph(K, HASHTABLE_SIZE, N_HT)
    kh2.consume_seqfile(filename2)

    for record in screed.open(paths):
        n = 0
        n_present = 0

        path = record.sequence
        n = len(path) - K + 1
        for i in range(n):
            if kh1.get(path[i:i + K]):
                n_present += 1

        if n_present / float(n) >= THRESHOLD:
            present1 = True
        else:
            present1 = False

        n = 0
        n_present = 0

        path = record.sequence
        n = len(path) - K + 1
        for i in range(n):
            if kh2.get(path[i:i + K]):
                n_present += 1

        if n_present / float(n) >= THRESHOLD:
            present2 = True
        else:
            present2 = False

        if present1 and not present2:
            print('>%s\n%s' % (record.name, record.sequence), file=uniq1)
        elif present2 and not present1:
            print('>%s\n%s' % (record.name, record.sequence), file=uniq2)


if __name__ == '__main__':
    main()
