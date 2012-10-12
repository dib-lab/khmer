import screed
import khmer
import sys

import math

hash_filename = sys.argv[1]
input_filename = sys.argv[2]
output_filename = sys.argv[3]

K=20 # word size
C=20 # use the C value from CBF post-diginorm

corrected = 0
miscorrected = 0
uncorrected = 0

outfp = open(output_filename, 'w')

ht = khmer.load_counting_hash(hash_filename)
aligner = khmer.new_readaligner(ht, 1, C)

for n, record in enumerate(screed.open(input_filename)):
   if n % 1000 == 0:
      print n

   seq = record.sequence
   name_tokens = record.name.split('_')
   seq_name = name_tokens[0]
   genome_seq = name_tokens[1]

   seq = seq.replace('-', '')
   genome_seq = genome_seq.replace('-', '')

   grXreAlign, reXgrAlign, score = aligner.align(seq)

   if len(reXgrAlign) > 0:
      graph_seq = grXreAlign.replace('-', '')
      if graph_seq in genome_seq or genome_seq in graph_seq:
         graph_len = len(graph_seq)
         gen_len = len(genome_seq)
         #if math.fabs(graph_len - gen_len) > 3:
         #   print 'difference > 3'
         #   print graph_seq
         #   print genome_seq
         #   print '---'
         corrected +=1 
      else:
         miscorrected += 1
      outfp.write('>%s\n%s\n' % (seq_name, graph_seq))
   else:
      uncorrected += 1
      outfp.write('>%s\n%s\n' % (seq_name, seq))

print 'corrected', corrected
print 'miscorrected', miscorrected
print 'uncorrected', uncorrected

outfp.close()
