#! /usr/bin/env python
import screed, sys

min_length = int(sys.argv[1])

for record in screed.fasta.fasta_iter(open(sys.argv[2])):
   if len(record['sequence']) >= min_length:
      print '>%s\n%s' % (record['name'], record['sequence'],)
