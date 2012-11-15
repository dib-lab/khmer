import screed
import khmer
import sys

import math

hash_filename = sys.argv[1]
input_filename = sys.argv[2]
output_filename = sys.argv[3]
max_error_region = int(sys.argv[4])

K=20 # word size
C=10 # 20

corrected = 0
uncorrected = 0

outfp = open(output_filename, 'w')

ht = khmer.load_counting_hash(hash_filename)
aligner = khmer.new_readaligner(ht, 1, C, max_error_region)

for n, record in enumerate(screed.open(input_filename)):
   if n % 1000 == 0:
      print n

   seq = record.sequence
   seq_name = record.name

   seq = seq.replace('N', 'A')

   grXreAlign, reXgrAlign = aligner.align(seq)

   if len(reXgrAlign) > 0:
      graph_seq = grXreAlign.replace('-', '')
      corrected += 1
      outfp.write('>%s\n%s\n' % (seq_name, graph_seq))
   else:
      uncorrected += 1
      outfp.write('>%s\n%s\n' % (seq_name, seq))


print 'corrected', corrected
print 'uncorrected', uncorrected

outfp.close()
