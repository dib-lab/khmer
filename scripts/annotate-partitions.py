#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Annotate sequences with partition numbers.

% python scripts/annotate-partitions.py <pmap_file> <file1> [ <file2> ... ]

Partition-annotated sequences will be in <fileN>.part.

Use '-h' for parameter help.
"""

import os
import argparse

import khmer
from khmer.file_api import check_file_status, check_space

DEFAULT_K = 32


def main():
    parser = argparse.ArgumentParser(
        description="Annotate seqs with partitions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--ksize', '-k', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('graphbase')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    ksize = args.ksize
    filenames = args.input_filenames
    htable = khmer.new_hashbits(ksize, 1, 1)

    partitionmap_file = args.graphbase + '.pmap.merged'

    check_file_status(partitionmap_file)
    for _ in filenames:
        check_file_status(_)

    check_space(filenames)

    print 'loading partition map from:', partitionmap_file
    htable.load_partitionmap(partitionmap_file)

    for infile in filenames:
        print 'outputting partitions for', infile
        outfile = os.path.basename(infile) + '.part'
        part_count = htable.output_partitions(infile, outfile)
        print 'output %d partitions for %s' % (part_count, infile)
        print 'partitions are in', outfile

if __name__ == '__main__':
    main()
