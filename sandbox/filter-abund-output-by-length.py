#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Trim sequences at k-mers of the given abundance, based on the given counting
hash table.  Output sequences will be placed in 'infile.abundfilt'.

% python scripts/filter-abund.py <counting.kh> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import screed.fasta
import os
import khmer
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader

from khmer.counting_args import build_counting_multifile_args

###

DEFAULT_CUTOFF = 2


class OutputByLength(object):
    def __init__(self, base):
        self.base = base
        self.fp_dict = {}

    def write(self, s):
        loc = s.find('\n')
        loc2 = s.find('\n', loc + 1)

        assert loc > -1
        assert loc2 > -1

        length = loc2 - loc
        assert length > 0

        fp_dict = self.fp_dict
        if length not in fp_dict:
            fp_dict[length] = open('%s.%03d' % (self.base, 1000 - length), 'w')

        fp_dict[length].write(s)


def main():
    parser = build_counting_multifile_args()
    parser.add_argument('--cutoff', '-C', dest='cutoff',
                        default=DEFAULT_CUTOFF, type=int,
                        help="Trim at k-mers below this abundance.")
    args = parser.parse_args()

    counting_ht = args.input_table
    infiles = args.input_filenames

    print 'file with ht: %s' % counting_ht

    print 'loading hashtable'
    ht = khmer.load_counting_hash(counting_ht)
    K = ht.ksize()

    print "K:", K

    ### the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = ht.trim_on_abundance(seq, args.cutoff)

        if trim_at >= K:
            return name, trim_seq

        return None, None

    ### the filtering loop
    for infile in infiles:
        print 'filtering', infile
        outfile = os.path.basename(infile) + '.abundfilt'

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), OutputByLength(outfile))

        print 'output in', outfile

if __name__ == '__main__':
    main()
