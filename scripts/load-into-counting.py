import sys, screed, os
import khmer

K = 32
HT_SIZE=2e9
N_HT=4

###

def main():
    base = sys.argv[1]
    filenames = sys.argv[2:]
    
    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    for n, filename in enumerate(filenames):
       print 'consuming input', filename
       ht.consume_fasta(filename)

       if n > 0 and n % 10 == 0:
          ht.save(base + '.ht')
          open(base + '.ht.info', 'w').write('through %s' % filename)

    ht.save(base + '.ht')

if __name__ == '__main__':
    main()
