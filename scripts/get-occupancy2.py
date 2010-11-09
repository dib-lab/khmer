import khmer, sys
import os.path

K=32
HASHTABLE_SIZE=int(1e9)

ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 1)

def callback(name, total_reads, n_consumed):
   print name, total_reads, n_consumed, ht.n_occupied()

def main(filename):
    basename = os.path.basename(filename)

    # populate the hash table and tag set
    ht.consume_fasta_and_tag(filename, callback)

if __name__ == '__main__':
    main(sys.argv[1])
