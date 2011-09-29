import sys, khmer
import os.path
import screed

K=32

readsfile = sys.argv[1]
contigfile = sys.argv[2]
outfile = os.path.basename(readsfile) + '.pathpart'
if len(sys.argv) == 4:
    outfile = sys.argv[3]

ht = khmer.new_hashbits(K, 1, 1)

# tag every 10th k-mer in the contigs
ht._set_tag_density(10)

ht.consume_fasta_and_tag(readsfile)

for n, record in enumerate(screed.open(contigfile)):
    if n % 10000 == 0:
        print '... joining', n
    ht.join_partitions_by_path(record.sequence)

#for n, record in enumerate(screed.open(readsfile)):
#    if n % 10000 == 0:
#        print '... joining x 2', n
#    ht.join_partitions_by_path(record.sequence)

print ht.output_partitions(readsfile, outfile)
