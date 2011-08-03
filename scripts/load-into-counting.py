# /usr/bin/env python
"""
Build a counting Bloom filter from the given sequences, save in <htname>.

% python scripts/load-into-counting.py <htname> <data1> [ <data2> <...> ]

Parameters to adjust: K, HT_SIZE.  HT_SIZE should be set to about 1/4 of the
available system memory.
"""

import sys, screed, os
import khmer

K = 32
HT_SIZE=2e7
N_HT=4

###

def main():
    base = sys.argv[1]
    filenames = sys.argv[2:]
    
    print 'making hashtable'
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)
    ht.set_use_bigcount(True)

    for n, filename in enumerate(filenames):
       print 'consuming input', filename
       ht.consume_fasta(filename)

       if n > 0 and n % 10 == 0:
           print 'mid-save', base
           ht.save(base)
           open(base + '.info', 'w').write('through %s' % filename)

    print 'saving', base
    ht.save(base)
    open(base + '.info', 'w').write('through end: %s' % filename)

if __name__ == '__main__':
    main()
