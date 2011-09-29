import sys, screed, os
import khmer

K = 32
N_HT=4

###

def main():
    ht_filename = sys.argv[1]
    contig_filename = sys.argv[2]

    print>>sys.stderr, 'loading ht from', ht_filename
    ht = khmer.new_counting_hash(K, 1, N_HT)
    ht.load(ht_filename)

    for record in screed.open(contig_filename):
       seq = record.sequence.upper()
       if 'N' in seq:
          seq = seq.replace('N', 'G')

       if K <= len(seq):
          a, b, c = ht.get_median_count(seq)
          print record.name, a, b, c, len(seq)
       else:
          print>>sys.stderr, 'skipping very short sequence', record.name, 'of length', len(seq)

if __name__ == '__main__':
    main()
