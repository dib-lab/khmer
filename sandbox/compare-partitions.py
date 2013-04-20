#! /usr/bin/env python
## only works for small files...

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

(filename1, filename2) = sys.argv[1:]

p1 = {}
s1 = {}
for name, pid, _, _ in read_partition_file(open(filename1)):
    name = name.split('\t')[0]
    x = p1.get(pid, set())
    x.add(name)
    p1[pid] = x

    s1[name] = pid

p2 = {}
s2 = {}
for name, pid, _, _ in read_partition_file(open(filename2)):
    name = name.split('\t')[0]
    x = p2.get(pid, set())
    x.add(name)
    p2[pid] = x

    s2[name] = pid

found = set()
for name in s1:
    pid = s1[name]
    pid2 = s2[name]

    x1 = p1[pid]
    x2 = p2[pid2]

    if x1 != x2 and pid not in found:
        print pid, pid2, len(x1), len(x2)
        found.add(pid)
