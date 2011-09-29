import sys
import khmer, screed

K = 32

ht = khmer.load_counting_hash(sys.argv[1])

for record in screed.open(sys.argv[2]):
   for pos in range(len(record.sequence) - K + 1):
      if ht.get(record.sequence[pos:pos+K]) > 40000:
         print record.sequence[pos:pos+K]
         sys.exit(0)
