import sys, khmer
import os.path
import screed

K=32

readsfile = sys.argv[1]
contigfile = sys.argv[2]
outfile = os.path.basename(readsfile) + '.path2part'
if len(sys.argv) == 4:
    outfile = sys.argv[3]

# create a hashbits data structure
ht = khmer.new_hashbits(K, 1, 1)

# tag every k-mer in the contigs
ht._set_tag_density(1)

# load contigs
ht.consume_fasta_and_tag(contigfile)

# iterate over the reads and join all the tags by overlapping read
for n, record in enumerate(screed.open(readsfile)):
    if n % 10000 == 0:
        print '... joining', n
    if 'N' in record.sequence:
        continue

    ht.join_partitions_by_path(record.sequence)

# output partitions
print ht.output_partitions(readsfile, outfile)
