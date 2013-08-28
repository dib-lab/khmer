#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import khmer
from screed.fasta import fasta_iter

K = 32

###


def get_partition(record):
    pid = record['name'].rsplit('\t', 1)[1]
    return int(pid)

###

ht = khmer.new_hashbits(K, 1, 1)

ht.consume_partitioned_fasta(sys.argv[1])
before = ht.count_partitions()

last_name = None
last_record = None
for n, record in enumerate(
        fasta_iter(open(sys.argv[1]), parse_description=False)):
    if n % 10000 == 0:
        print '...', n

    name = record['name'].split()[0]
    name = name.split('/', 1)[0]

    if name == last_name:
        if 1:
            pid_1 = ht.get_partition_id(last_record['sequence'][:K])
            pid_2 = ht.get_partition_id(record['sequence'][:K])

            ht.join_partitions(pid_1, pid_2)
        else:                           # TEST
            pid_1 = get_partition(last_record)
            pid_2 = get_partition(record)
            assert pid_1 == pid_2, (last_record, record, pid_1, pid_2)

    last_name = name
    last_record = record

ht.output_partitions(sys.argv[1], sys.argv[1] + '.paired')
print 'before:', before
after = ht.count_partitions()
print 'after:', after

n_combined = before[0] - after[0]
print 'combined:', n_combined

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
