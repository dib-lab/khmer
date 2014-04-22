#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys

MAX_SIZE = 50000


def read_partition_file(fp):
    for n, line in enumerate(fp):
        if n % 2 == 0:
            surrendered = False
            name, partition_id, partition_readcount = line.split('\t')

            partition_readcount = int(partition_readcount)
            if '*' in partition_id:
                partition_id = int(partition_id[:-1])
                surrendered = True
            else:
                partition_id = int(partition_id)
        else:
            sequence = line.strip()

            yield name, partition_id, partition_readcount, surrendered, \
                sequence

if __name__ == '__main__':
    filename = sys.argv[1]
    prefix = sys.argv[2]

    partition_sizes = {}

    # first, read in all the cluster sizes

    fp = open(filename)
    for n, x in enumerate(read_partition_file(fp)):
        if n % 100000 == 0:
            print '...', n

        name, partition_id, readcount, surrendered, seq = x
        if not surrendered:
            partition_sizes[partition_id] = readcount

    # sort by # of reads in each cluster
    divvy = sorted(partition_sizes.items(), key=lambda y: y[1])

    # divvy up into different groups, based on having MAX_SIZE sequences
    # in each group.
    total = 0
    group = set()
    group_n = 0
    group_d = {}
    for partition_id, n_reads in divvy:
        group.add(partition_id)
        total += n_reads

        if total > MAX_SIZE:
            for partition_id in group:
                group_d[partition_id] = group_n

            group_n += 1
            group = set()
            total = 0

    print '%d groups' % group_n

    # open a bunch of output files for the different groups
    group_fps = {}
    for n in range(group_n):
        fp = open('%s.group%d.fa' % (prefix, n), 'w')
        group_fps[n] = fp

    surrendered_fp = open('%s.surrender.fa' % prefix, 'w')

    # write 'em all out!
    fp = open(filename)
    for n, x in enumerate(read_partition_file(fp)):
        if n % 100000 == 0:
            print '...x2', n

        name, partition_id, readcount, surrendered, seq = x
        if surrendered:
            outfp = surrendered_fp
        else:
            group_n = group_d[partition_id]
            outfp = group_fps[group_n]

        outfp.write('>%s %s %s\n%s\n' % (name, partition_id, readcount, seq))
