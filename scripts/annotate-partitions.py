#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
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

#  Import fileapi from sandbox - temporary arrangement
current_file_path = os.path.realpath(__file__)
current_folder = os.path.dirname(current_file_path)
parent_folder = os.path.dirname(current_folder)
sandbox_folder = os.path.join(parent_folder, 'sandbox')
sys.path.append(sandbox_folder)

import fileApi

DEFAULT_K = 32


def main():
    parser = argparse.ArgumentParser(
        description="Annotate seqs with partitions.")

    parser.add_argument('--ksize', '-k', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('graphbase')
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    K = args.ksize
    ht = khmer.new_hashbits(K, 1, 1)

    partitionmap_file = args.graphbase + '.pmap.merged'

    print 'loading partition map from:', partitionmap_file
    ht.load_partitionmap(partitionmap_file)

    # Check space availability
    freeSpace = fileApi.check_space(args.input_filenames)
    if freeSpace != 0:
        print >>sys.stderr, 'ERROR: Not enough free space on disk, \
        need at least %s more,' % str(freeSpace)
        sys.exit(-1)
    
    for infile in args.input_filenames:
        print 'outputting partitions for', infile
        outfile = os.path.basename(infile) + '.part'
        n = ht.output_partitions(infile, outfile)
        print 'output %d partitions for %s' % (n, infile)
        print 'partitions are in', outfile

if __name__ == '__main__':
    main()
