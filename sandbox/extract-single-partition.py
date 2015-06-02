#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
from __future__ import print_function
import sys
from screed.fasta import fasta_iter


def read_partition_file(fp):
    for n, record in enumerate(fasta_iter(fp, parse_description=False)):
        name = record['name']
        name, partition_id = name.rsplit('\t', 1)
        yield n, name, int(partition_id), record['sequence']


def main():
    select_pid = int(sys.argv[2])
    count = 0
    for n, name, pid, seq in read_partition_file(open(sys.argv[1])):
        if pid == select_pid:
            print('>%s\t%d\n%s' % (name, pid, seq))
            count += 1

        if n % 10000 == 0:
            sys.stderr.write('...%d\n' % (n,))

    sys.stderr.write('found %d total in partition %d\n' % (count, pid))


if __name__ == '__main__':
    main()
