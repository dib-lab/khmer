#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Extract partitioned sequences into files grouped by partition size.

% python scripts/extract-partitions.py <base> <file1.part> [ <file2.part> ... ]

Grouped sequences will be <base>.groupN.fa files.

Use '-h' for parameter help.

@CTB note that if threshold is != 1, those sequences will not be output
by output_unassigned...
"""

import sys
import screed
import argparse
from khmer.file import check_file_status, check_space

DEFAULT_MAX_SIZE = int(1e6)
DEFAULT_THRESHOLD = 5


def read_partition_file(filename):
    for record_index, record in enumerate(screed.open
                                          (filename, parse_description=False)):
        name, partition_id = record.name.rsplit('\t', 1)
        yield record_index, record, int(partition_id)


def output_single(read):
    if hasattr(read, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.accuracy)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)

#


def main():
    parser = argparse.ArgumentParser(
        description="Extract partitioned seqs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

    max_size = args.max_size
    threshold = args.min_part_size
    output_groups = args.output_groups
    output_unassigned = args.output_unass

    prefix = args.prefix
    distfilename = prefix + '.dist'

    n_unassigned = 0

    for infile in args.part_filenames:
        check_file_status(infile)

    check_space(args.part_filenames)

    print '---'
    print 'reading partitioned files:', repr(args.part_filenames)
    if output_groups:
        print 'outputting to files named "%s.groupN.fa"' % prefix
        print 'min reads to keep a partition:', threshold
        print 'max size of a group file:', max_size
    else:
        print 'NOT outputting groups! Beware!'

    if output_unassigned:
        print 'outputting unassigned reads to "%s.unassigned.fa"' % prefix

    print 'partition size distribution will go to %s' % distfilename
    print '---'

    #

    suffix = 'fa'
    is_fastq = False

    for index, read, pid in read_partition_file(args.part_filenames[0]):
        if hasattr(read, 'accuracy'):
            suffix = 'fq'
            is_fastq = True
        break

    for filename in args.part_filenames:
        for index, read, pid in read_partition_file(filename):
            if is_fastq:
                assert hasattr(read, 'accuracy'), \
                    "all input files must be FASTQ if the first one is"
            else:
                assert not hasattr(read, 'accuracy'), \
                    "all input files must be FASTA if the first one is"

            break

    if output_unassigned:
        unassigned_fp = open('%s.unassigned.%s' % (prefix, suffix), 'w')

    count = {}
    for filename in args.part_filenames:
        for index, read, pid in read_partition_file(filename):
            if index % 100000 == 0:
                print '...', index

            count[pid] = count.get(pid, 0) + 1

            if pid == 0:
                n_unassigned += 1
                if output_unassigned:
                    print >>unassigned_fp, output_single(read)

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
    for counter, index in sorted(dist.items()):
        total += index
        wtotal += counter * index
        distfp.write('%d %d %d %d\n' % (counter, index, total, wtotal))
    distfp.close()

    if not output_groups:
        sys.exit(0)

    # sort groups by size
    divvy = sorted(count.items(), key=lambda y: y[1])
    divvy = filter(lambda y: y[1] > threshold, divvy)

    # divvy up into different groups, based on having max_size sequences
    # in each group.
    total = 0
    group = set()
    group_n = 0
    group_d = {}
    for partition_id, n_reads in divvy:
        group.add(partition_id)
        total += n_reads

        if total > max_size:
            for partition_id in group:
                group_d[partition_id] = group_n
                # print 'group_d', partition_id, group_n

            group_n += 1
            group = set()
            total = 0

    if group:
        for partition_id in group:
            group_d[partition_id] = group_n
            # print 'group_d', partition_id, group_n
        group_n += 1

    print '%d groups' % group_n
    if group_n == 0:
        print 'nothing to output; exiting!'
        return

    # open a bunch of output files for the different groups
    group_fps = {}
    for _ in range(group_n):
        group_fp = open('%s.group%04d.%s' % (prefix, _, suffix), 'w')
        group_fps[_] = group_fp

    # write 'em all out!

    total_seqs = 0
    part_seqs = 0
    toosmall_parts = 0
    for filename in args.part_filenames:
        for index, read, partition_id in read_partition_file(filename):
            total_seqs += 1
            if index % 100000 == 0:
                print '...x2', index

            if partition_id == 0:
                continue

            try:
                group_n = group_d[partition_id]
            except KeyError:
                assert count[partition_id] <= threshold
                toosmall_parts += 1
                continue

            outfp = group_fps[group_n]

            outfp.write(output_single(read))
            part_seqs += 1

    print '---'
    print 'Of %d total seqs,' % total_seqs
    print 'extracted %d partitioned seqs into group files,' % part_seqs
    print 'discarded %d sequences from small partitions (see -m),' % \
        toosmall_parts
    print 'and found %d unpartitioned sequences (see -U).' % n_unassigned
    print ''
    print 'Created %d group files named %s.groupXXXX.%s' % (len(group_fps),
                                                            prefix,
                                                            suffix)

if __name__ == '__main__':
    main()
