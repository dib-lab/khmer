#! /usr/bin/env python
"""
Extract partitioned sequences into files grouped by partition size.

% python scripts/extract-partitions.py <base> <file1.part> [ <file2.part> ... ]

Grouped sequences will be <base>.groupN.fa files.

Use '-h' for parameter help.

@CTB note that if THRESHOLD is != 1, those sequences will not be output
by output_unassigned...
"""

import sys
import os.path
import screed
import argparse

DEFAULT_MAX_SIZE=int(1e5)
DEFAULT_THRESHOLD=5

def read_partition_file(filename):
    for n, record in enumerate(screed.open(filename)):
        name = record.name
        name, partition_id = name.rsplit('\t', 1)
        yield n, name, int(partition_id), record.sequence

###

def main():
    parser = argparse.ArgumentParser(description="Extract partitioned seqs.")

    parser.add_argument('prefix')
    parser.add_argument('part_filenames', nargs='+')
    parser.add_argument('--max-size', '-X', dest='max_size',
                        default=DEFAULT_MAX_SIZE, type=int,
                        help='Max group size (n sequences)')
    parser.add_argument('--min-partition-size', '-m', dest='min_part_size',
                        default=DEFAULT_THRESHOLD, type=int,
                        help='Minimum partition size worth keeping')
    parser.add_argument('--no-output-groups', '-n', dest='output_groups',
                        default=True, action='store_false',
                        help='Do not actually output groups files.')
    parser.add_argument('--output-unassigned', '-U', dest='output_unass',
                        default=False, action='store_true',
                        help='Output unassigned sequences, too')

    args = parser.parse_args()

    MAX_SIZE = args.max_size
    THRESHOLD = args.min_part_size
    output_groups = args.output_groups
    output_unassigned = args.output_unass

    prefix = args.prefix
    distfilename = prefix + '.dist'

    print '---'
    print 'reading partitioned files:', repr(args.part_filenames)
    if output_groups:
        print 'outputting to files named "%s.groupN.fa"' % prefix
        print 'min reads to keep a partition:', THRESHOLD
        print 'max size of a group file:', MAX_SIZE
    else:
        print 'NOT outputting groups! Beware!'

    if output_unassigned:
        print 'outputting unassigned reads to "%s.unassigned.fa"' % prefix

    print 'partition size distribution will go to %s' % distfilename
    print '---'

    ###

    if output_unassigned:
        unassigned_fp = open('%s.unassigned.fa' % prefix, 'w')

    count = {}
    for filename in args.part_filenames:
        for n, name, pid, seq in read_partition_file(filename):
            if n % 100000 == 0:
                print '...', n

            count[pid] = count.get(pid, 0) + 1

            if pid == 0 and output_unassigned:
                print >>unassigned_fp, '>%s\n%s' % (name, seq)

    if output_unassigned:
        unassigned_fp.close()

    if 0 in count:                          # eliminate unpartitioned sequences
       del count[0]

    # develop histogram of partition sizes
    dist = {}
    for pid, size in count.items():
        dist[size] = dist.get(size, 0) + 1

    # output histogram
    distfp = open(distfilename, 'w')

    total = 0
    wtotal = 0
    for c, n in sorted(dist.items()):
        total += n
        wtotal += c*n
        distfp.write('%d %d %d %d\n' % (c, n, total, wtotal))
    distfp.close()

    if not output_groups:
        sys.exit(0)

    # sort groups by size
    divvy = sorted(count.items(), key=lambda y:y[1])
    divvy = filter(lambda y:y[1] > THRESHOLD, divvy)

    ## divvy up into different groups, based on having MAX_SIZE sequences
    ## in each group.
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
                #print 'group_d', partition_id, group_n

            group_n += 1
            group = set()
            total = 0

    if group:
        for partition_id in group:
            group_d[partition_id] = group_n
            #print 'group_d', partition_id, group_n
        group_n += 1


    print '%d groups' % group_n
    if group_n == 0:
        print 'nothing to output; exiting!'
        return

    ## open a bunch of output files for the different groups
    group_fps = {}
    for n in range(group_n):
        fp = open('%s.group%04d.fa' % (prefix, n), 'w')
        group_fps[n] = fp

    ## write 'em all out!

    for filename in args.part_filenames:
        for n, name, partition_id, seq in read_partition_file(filename):
            if n % 100000 == 0:
                print '...x2', n

            if partition_id == 0:
                continue

            try:
                group_n = group_d[partition_id]
            except KeyError:
                assert count[partition_id] <= THRESHOLD
                continue

            outfp = group_fps[group_n]

            outfp.write('>%s\t%s\n%s\n' % (name, partition_id, seq))

if __name__ == '__main__':
    main()
