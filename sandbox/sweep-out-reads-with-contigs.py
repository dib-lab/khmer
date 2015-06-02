#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import sys
import khmer
import os.path
import screed

K = 20


def main():
    readsfile = sys.argv[1]
    contigfile = sys.argv[2]
    outfile = os.path.basename(readsfile) + '.sweep'
    if len(sys.argv) == 4:
        outfile = sys.argv[3]

    # create a hashbits data structure
    ht = khmer.new_hashbits(K, 1, 1)

    # tag every k-mer in the contigs
    ht._set_tag_density(0)

    # load contigs, connect into N partitions
    print('loading contigs from', contigfile)
    ht.consume_fasta_and_tag(contigfile)
    subset = ht.do_subset_partition(0, 0)
    ht.merge_subset(subset)

    print('outputting contig-partitioned reads to', outfile)
    ht.output_partitions(readsfile, outfile, True)


if __name__ == '__main__':
    main()
