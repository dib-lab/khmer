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
"""
Trim sequences at k-mers of the given abundance, based on the given counting
hash table.  Output sequences will be placed in 'infile.abundfilt'.

% python sandbox/split-sequences-by-length.py <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

###

DEFAULT_CUTOFF = 2


class OutputByLength(object):

    def __init__(self, base):
        self.base = base
        self.fp_dict = {}

    def save(self, name, sequence):
        length = len(sequence)

        fp_dict = self.fp_dict
        if length not in fp_dict:
            fp_dict[length] = open('%s.%03d' % (self.base, 1000 - length), 'w')

        fp_dict[length].write('>%s\n%s\n' % (name, sequence))


def main():
    base = sys.argv[1]
    filenames = sys.argv[2:]

    out = OutputByLength(base)

    n = 0
    for filename in filenames:
        print('opening')
        for record in screed.open(filename):
            out.save(record.name, record.sequence)
            n += 1
            if n % 10000 == 0:
                print('...', n)

if __name__ == '__main__':
    main()
