#! /usr/bin/env python
import sys
from screed.fasta import fasta_iter

MAX_SIZE=5000

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


threshold = int(sys.argv[3])

prefix = sys.argv[2]
distfp = open(sys.argv[4], 'w')
count = {}
filename = sys.argv[1]

for n, record in enumerate(fasta_iter(open(sys.argv[1]))):
    if n % 10000 == 0:
        print '...', n

    partition = int(record['name'].rsplit('\t', 1)[1].rstrip('*'))
    count[partition] = count.get(partition, 0) + 1

# develop histogram of partition sizes
dist = {}
for n, record in enumerate(fasta_iter(open(sys.argv[1]))):
    if n % 10000 == 0:
        print '...x2', n

    partition = int(record['name'].rsplit('\t', 1)[1].rstrip('*'))
    if partition not in count:
        continue
    
    c = count[partition]
    if partition == 0:
        c = 0

    dist[c] = dist.get(c, 0) + 1

#    if c >= threshold:
#        outfp.write('>%s\n%s\n' % (record['name'], record['sequence']))

# output histogram
total = 0
for c, n in sorted(dist.items()):
    if c:
        n /= c
    total += n
    distfp.write('%d %d %d\n' % (c, n, total))

# separate
if 0 in count:
   del count[0]

divvy = sorted(count.items(), key=lambda y:y[1])

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
            print 'group_d', partition_id, group_n

        group_n += 1
        group = set()
        total = 0

if group:
    for partition_id in group:
        group_d[partition_id] = group_n
        print 'group_d', partition_id, group_n
    group_n += 1


print '%d groups' % group_n

## open a bunch of output files for the different groups
group_fps = {}
for n in range(group_n):
    fp = open('%s.group%d.fa' % (prefix, n), 'w')
    group_fps[n] = fp

surrendered_fp = open('%s.surrender.fa' % prefix, 'w')

## write 'em all out!
fp = open(filename)
for n, x in enumerate(read_partition_file(fp)):
    if n % 100000 == 0:
        print '...x2', n

    name, partition_id, surrendered, seq = x
    if partition_id == 0:
        continue
    
    if surrendered:
        surrender_ch = '*'
        outfp = surrendered_fp
    else:
        surrender_ch = ' '
        group_n = group_d[partition_id]
        outfp = group_fps[group_n]

    outfp.write('>%s\t%s%s\n%s\n' % (name, partition_id, surrender_ch, seq))
