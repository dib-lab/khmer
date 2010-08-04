import sys
import glob
import khmer

HASHTABLE_SIZE = 4**12+1                # in bytes; default, 17gb

files = glob.glob(sys.argv[1])
K = int(sys.argv[2])

total_reads = len(files)*[0]
n_consumed = len(files)*[0]
n_seq_kept = len(files)*[0]

print 'allocating hashtable: ', HASHTABLE_SIZE
ht = khmer.new_hashtable(K, HASHTABLE_SIZE)

for i, infile in enumerate(files):
   (total_reads[i], n_consumed[i]) = ht.consume_fasta(infile)

for i, infile in enumerate(files):
   ht.output_fasta_kmer_pos_freq(infile, infile + "." + str(K)  + ".freq")
