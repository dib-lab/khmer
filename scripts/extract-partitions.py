#! /usr/bin/env python
import sys
from screed.fasta import fasta_iter

MAX_SIZE=50000
THRESHOLD=1

def read_partition_file(fp):
    for n, record in enumerate(fasta_iter(fp, parse_description=False)):
        name = record['name']
        name, partition_id = name.rsplit('\t', 1)
        yield n, name, partition_id, record['sequence']

###

(filename, prefix) = sys.argv[1:]

distfilename = prefix + '.dist'
distfp = open(distfilename, 'w')

print '---'
print 'reading partitioned file:', filename
print 'outputting to files named "%s.groupN.fa"' % prefix
print 'partition size distribution will go to %s' % distfilename
print '---'

###

count = {}

###

for n, name, pid, seq in read_partition_file(open(filename)):
    if n % 100000 == 0:
        print '...', n

    count[pid] = count.get(pid, 0) + 1

if 0 in count:                          # eliminate unpartitioned sequences
   del count[0]
   
# develop histogram of partition sizes
dist = {}
for pid, size in count.items():
    dist[size] = dist.get(size, 0) + 1
        
# output histogram
total = 0
for c, n in sorted(dist.items()):
    total += n
    distfp.write('%d %d %d\n' % (c, n, total))
distfp.close()

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

## open a bunch of output files for the different groups
group_fps = {}
for n in range(group_n):
    fp = open('%s.group%d.fa' % (prefix, n), 'w')
    group_fps[n] = fp

## write 'em all out!
    
fp = open(filename)
for n, name, partition_id, seq in read_partition_file(fp):
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
