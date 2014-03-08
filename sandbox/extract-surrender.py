#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import sys
from screed.fasta import fasta_iter


def read_partition_file(fp):
    for n, line in enumerate(fp):
        if n % 2 == 0:
            surrendered = False
            name, partition_id = line[1:].strip().rsplit('\t', 1)

            if '*' in partition_id:
                partition_id = int(partition_id[:-1])
                surrendered = True
            else:
                partition_id = int(partition_id)
        else:
            sequence = line.strip()

            yield name, partition_id, surrendered, sequence


def main(filename):
    outfp = open(filename + '.surrender', 'w')
    for n, x in enumerate(read_partition_file(open(filename))):
        name, part_id, surrendered, sequence = x
        if n % 10000 == 0:
            print '...', n

        if surrendered:
            outfp.write('>%s\n%s\n' % (name, sequence))

if __name__ == '__main__':
    filename = sys.argv[1]
    main(filename)
