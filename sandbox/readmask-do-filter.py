#! /usr/bin/env python
import sys
import khmer

sequence_file = sys.argv[1]
readmask_file = sys.argv[2]
output_file = sys.argv[3]

r = khmer.new_readmask(0)
r.load(readmask_file)
r.filter_fasta_file(sequence_file, output_file)
