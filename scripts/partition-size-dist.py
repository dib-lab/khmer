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


distfp = open(sys.argv[2], 'w')
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
