import sys
import khmer

files = sys.argv[2:]

total_reads = len(files) * [0]
n_consumed = len(files) * [0]
n_seq_kept = len(files) * [0]

print 'loading ht'
ht = khmer.new_counting_hash(1, 1, 1)

ht.load(sys.argv[1])

for i, infile in enumerate(files):
    ht.output_fasta_kmer_pos_freq(infile, infile + ".freq")
