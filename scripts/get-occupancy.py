import khmer
import sys
import screed

K = 32
#HASHTABLE_SIZE=int(256e9)
HASHTABLE_SIZE = 1000000000

print 'creating ht with size ' + str(HASHTABLE_SIZE)
ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

read_count = 0

fd = open(sys.argv[1], "w")

for filename in sys.argv[2:]:
   print "processing file: " + filename + " reads processed: " + str(read_count)
 
   for n, record in enumerate(screed.fastq.fastq_iter(open(filename))):
      read_count += 1
      ht.consume(record['sequence'])

      if read_count % 1000000 == 0:
         fd.write(str(read_count) + " " + str(ht.n_occupied()) + "\n")
         fd.flush()

fd.write(str(read_count) + " " + str(ht.n_occupied()) + "\n")
fd.close()

