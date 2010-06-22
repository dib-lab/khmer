import sys
import glob
import khmer

K = 17
THRESHOLD = 5
MINLENGTH = 50

infile = sys.argv[1]
#outfile = sys.argv[2]

hasher = khmer.new_hashtable(K, 16000000011)
hasher.consume_fasta(infile)

