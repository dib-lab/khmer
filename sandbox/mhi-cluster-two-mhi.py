#! /usr/bin/env python
"compare two mhi"
from __future__ import print_function
import khmer
import screed
import argparse
import os.path
try:
    import numpy
except ImportError:
    pass

KSIZE=32
COMBINED_MH_SIZE=1000

def load_and_tag(ct, filename):
    print('reading and tagging sequences')
    for record in screed.open(filename):
        print('.', record.name)
        ct.consume_and_tag(record.sequence)
    print('...done loading sequences!')

def filter_combined(combined, min_tagcount=100):
    return [ c for c in combined if len(c.get_tags()) >= min_tagcount ]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile1')
    parser.add_argument('infile2')
    parser.add_argument('prefix')
    args = parser.parse_args()
    
    infile1 = args.infile1
    infile2 = args.infile2
    
    print('loading nbhd minhashes 1...')
    nbhd_mh1 = khmer.load_neighborhood_minhash(infile1)
    print('...done!')
    total_tags = len(nbhd_mh1.get_all_tags())
    
    print('loading nbhd minhashes 2...')
    nbhd_mh2 = khmer.load_neighborhood_minhash(infile2)
    print('...done!')

    print('building ~chromosome level minhashes 1')
    combined1 = nbhd_mh1.build_combined_minhashes(COMBINED_MH_SIZE)
    combined1 = filter_combined(combined1)

    tags_in_combined1 = sum([ len(c.get_tags()) for c in combined1 ])
    print('xxx', total_tags, tags_in_combined1)

    print('building ~chromosome level minhashes 2')
    combined2 = nbhd_mh2.build_combined_minhashes(COMBINED_MH_SIZE)

    matched1 = set()
    matched2 = set()
    for i, x in enumerate(combined1):
        for j, y in enumerate(combined2):
            d = x.get_minhash().compare(y.get_minhash())
            if d > 0.05:
                matched1.add(x)
                matched2.add(y)

    print('in common: %d, %d' % (len(matched1), len(matched2)))

    compare = list(matched1)
    compare_names = ["s1.%d" % i for i in range(len(matched1)) ]
    compare += list(matched2)
    compare_names += ["s2.%d" % i for i in range(len(matched2)) ]
    
    D = numpy.zeros([len(compare), len(compare)])
    numpy.set_printoptions(precision=3, suppress=True)
    for i, x in enumerate(compare):
        for j, y in enumerate(compare):
            d = x.get_minhash().compare(y.get_minhash())
            if d > 0.05:
                matched1.add(x)
                matched2.add(y)
                D[i, j] = d

    print(D)
    fp = open(args.prefix, 'wb')
    numpy.save(fp, D)
    fp.close()

    fp = open('%s.labels.txt' % args.prefix, 'w')
    fp.write("\n".join(compare_names))
    fp.close()
    

if __name__ == '__main__':
    main()
    #nbhd_mh, combined = main()
