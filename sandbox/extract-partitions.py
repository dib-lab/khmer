#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import os.path
import screed
import argparse
from khmer.counting_args import build_construct_args

MAX_SIZE = int(1e6)
THRESHOLD = 1

output_groups = True
output_unassigned = False


def read_partition_file(filename):
    
    # For each record in input file, yield its index, name, partition ID and sequence
    # Input file format looks like:
    #
    # >name    partition_id
    # SEQUENCE
    #
    for n, record in enumerate(screed.open(filename)):
        name = record.name
        name, partition_id = name.rsplit('\t', 1)
        yield n, name, int(partition_id), record.sequence

###
parser = build_construct_args("Extract partitions")
parser.add_argument('filename')
parser.add_argument('-p','--prefix',dest='prefix',help='Prefix to group names')
args=parser.parse_args()
filename = args.filename
# If prefix is given as second arg, pick it up.
# Else, use filename-sans-path as prefix
prefix=args.prefix if args.prefix else os.path.basename(filename) 



# Create dist file and open to write
distfilename = prefix + '.dist'
distfp = open(distfilename, 'w')

print '---'
print 'reading partitioned file:', filename

# Output reads after grouping them by partition annotation
if output_groups:
    print 'outputting to files named "%s.groupN.fa"' % prefix
    print 'min reads to keep a partition:', THRESHOLD
    print 'max size of a group file:', MAX_SIZE
else:
    print 'NOT outputting groups! Beware!'

# Output reads unassigned to any partition to a separate file
if output_unassigned:
    print 'outputting unassigned reads to "%s.unassigned.fa"' % prefix

# Output details on partition size distribution
print 'partition size distribution will go to %s' % distfilename
print '---'

###

count = {}

###

# Open unassigned-reads-file to write into
if output_unassigned:
    unassigned_fp = open('%s.unassigned.fa' % prefix, 'w')

# Get details of each sequence
for n, name, pid, seq in read_partition_file(filename):
    # Progress indicator
    if n % 100000 == 0:
        print '...', n

    # Count number of sequences in partition with ID=pid
    # If this is the first member of a partition, initialize
    # count to 1
    count[pid] = count.get(pid, 0) + 1

    # If current pid is 0 (no partition allocated), write
    # read to unassigned-reads-file
    if pid == 0 and output_unassigned:
        print >>unassigned_fp, '>%s\n%s' % (name, seq)

# Close unassigned-reads-file if open
if output_unassigned:
    unassigned_fp.close()

if 0 in count:                          # eliminate unpartitioned sequences
    del count[0]

# develop histogram of partition sizes
# operation equivalent to 
# COUNT(*) GROUP BY PARTITION_SIZE  
dist = {}
for pid, size in count.items():
    dist[size] = dist.get(size, 0) + 1

# output histogram
total = 0
wtotal = 0
# Iterate thru sorted group-size-count distribution and
# write group-size, group-size-count, cumul-sum, weighted-cumul-sum
for c, n in sorted(dist.items()):
    total += n
    wtotal += c * n
    distfp.write('%d %d %d %d\n' % (c, n, total, wtotal))
distfp.close()

if not output_groups:
    sys.exit(0)

# sort groups by size and filter out groups with size
# below THRESHOLD
divvy = sorted(count.items(), key=lambda y: y[1])
divvy = filter(lambda y: y[1] > THRESHOLD, divvy)

## divvy up into different groups, based on having MAX_SIZE sequences
## in each group.
total = 0
group = set()
group_n = 0
group_d = {}
# Iterating through pid:count(reads) in sorted+filtered group-data,
for partition_id, n_reads in divvy:
    # We add the pid to our "group" set;
    group.add(partition_id)
    # and the count(reads) to a total counter
    total += n_reads
    
    # If the total count of reads in groups tracked until current iteration
    # is more than a pre-determined max size 
    if total > MAX_SIZE:
        # package the current groups and open a new container for
        # the next set of groups. Each pack bears the ID of the
        # last partition of the group of partitions
        for partition_id in group:
            group_d[partition_id] = group_n
            # print 'group_d', partition_id, group_n

        group_n += 1
        group = set()
        total = 0
# Final cleanup - for the final batch of groups processed,
# because their total read count might not be >MAX_SIZE 
if group:
    for partition_id in group:
        group_d[partition_id] = group_n
        # print 'group_d', partition_id, group_n
    group_n += 1


print '%d groups' % group_n

## open a bunch of output files for the different groups
group_fps = {}
for n in range(group_n):
    fp = open('%s.group%04d.fa' % (prefix, n), 'w')
    group_fps[n] = fp

## write 'em all out!

# Get index, seqName, pid and seq from inFile and process as
# follows:
for n, name, partition_id, seq in read_partition_file(filename):
    # If 100,000 seqs have been processed, print a progress message
    if n % 100000 == 0:
        print '...x2', n

    # If no partition allocated to read, skip it
    if partition_id == 0:
        continue

    # Check if there exists a group of partition with the current
    # pid
    try:
        group_n = group_d[partition_id]
    except KeyError:
        # If not, ensure that count of seqs with the matching pid
        # is > THRESHOLD and skip
        assert count[partition_id] <= THRESHOLD
        continue

    # Pick the current output file
    outfp = group_fps[group_n]

    # Write to the output file
    outfp.write('>%s\t%s\n%s\n' % (name, partition_id, seq))
