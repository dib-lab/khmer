import string
import screed
import khmer
import sys

hash_filename = sys.argv[1]
input_filename = sys.argv[2]
#output_filename = sys.argv[3]

K=20

#outfp = open(output_filename, 'w')

ht = khmer.load_counting_hash(hash_filename)
aligner = khmer.new_readaligner(ht)

base_counter = 0
false_positives = 0
errors_found = 0
errors_fixed = 0

for n, record in enumerate(screed.open(input_filename)):
   seq = record.sequence
   base_counter += len(seq)

   for base in seq:
      if base in string.ascii_lowercase:
         errors_found += 1

   graphAlign, readAlign, score = aligner.align(seq)

   if len(readAlign) > 0:
      readMap = []
      count = 0
      for baseNo, base in enumerate(readAlign):
         if base != '-':
            readMap.append(count)
         count += 1

      for baseNo, base in enumerate(seq):
         corr_base_index = readMap[baseNo]

         is_error = False
         if base in string.ascii_lowercase:
            is_error = True
            base = base.upper()

            if is_error:
               if base != graphAlign[corr_base_index]:
                  errors_fixed += 1
            else:
               if base != graphAlign[corr_base_index]:
                  false_positives += 1

print 'errors_found: ', errors_found
print 'errors_fixed: ', errors_fixed
print 'false_positives: ', false_positives
print 'base_counter: ', base_counter

