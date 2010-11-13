import khmer
import sys
import screed
import os
import subprocess
import glob
import zlib

K = 32
HASHTABLE_SIZE=int(512e9)
#HASHTABLE_SIZE = 1000000000

print 'creating ht with size ' + str(HASHTABLE_SIZE)
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

read_count = 0

fd = open(sys.argv[1], "w")

files = glob.glob(sys.argv[2])

for filename in files:
   print "processing file: " + filename + " reads processed: " + str(read_count)
 
   for n, record in enumerate(screed.fastq.fastq_iter(gzip.open(filename))):
      read_count += 1
      ht.consume(record['sequence'])

      if read_count % 1000000 == 0:
         fd.write(str(read_count) + " " + str(ht.n_occupied()) + "\n")
         fd.flush()

fd.write(str(read_count) + " " + str(ht.n_occupied()) + "\n")
fd.close()

