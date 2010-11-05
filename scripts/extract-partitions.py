#! /usr/bin/env python
import sys
from screed.fasta import fasta_iter

MAX_SIZE=50000
THRESHOLD=1

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


(filename, prefix, distfilename) = sys.argv[1:]

distfp = open(distfilename, 'w')
surrendered_exist = False

count = {}

###

for n, (name, pid, surr, seq) in enumerate(read_partition_file(open(filename))):
    if n % 10000 == 0:
        print '...', n

    count[pid] = count.get(pid, 0) + 1
    if surr:
        surrendered_exist = True

# develop histogram of partition sizes
dist = {}
for n, (name, pid, surr, seq) in enumerate(read_partition_file(open(filename))):
    if n % 10000 == 0:
        print '...x2', n

    if pid not in count:
        continue
    
    c = count[pid]
    if pid == 0:
        c = 0

    dist[c] = dist.get(c, 0) + 1

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

if surrendered_exist:
    surrendered_fp = open('%s.surrender.fa' % prefix, 'w')

## write 'em all out!
fp = open(filename)
for n, x in enumerate(read_partition_file(fp)):
    if n % 100000 == 0:
        print '...x3', n

    name, partition_id, surrendered, seq = x
    if partition_id == 0:
        continue
    
    if surrendered:
        surrender_ch = '*'
        outfp = surrendered_fp
    else:
        surrender_ch = ' '
        try:
            group_n = group_d[partition_id]
        except KeyError:
            continue
        outfp = group_fps[group_n]

    outfp.write('>%s\t%s%s\n%s\n' % (name, partition_id, surrender_ch, seq))
