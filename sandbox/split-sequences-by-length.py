#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Trim sequences at k-mers of the given abundance, based on the given counting
hash table.  Output sequences will be placed in 'infile.abundfilt'.

% python scripts/filter-abund.py <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
from __future__ import print_function
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
