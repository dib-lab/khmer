#! /usr/bin/env python
"search mhi nbhd"
from __future__ import print_function
import khmer
import screed
from cPickle import dump
import argparse
import os.path
import numpy

import sys
sys.path.append('../sourmash')
try:
    import sourmash_lib, sourmash_signature
except ImportError:
    pass

KSIZE=32
COMBINED_MH_SIZE=500
COMBINE_THIS_MANY=10000

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('sigfile')
    args = parser.parse_args()
    
    infile = args.infile
    
    print('loading nbhd minhashes...')
    nbhd_mh = khmer._minhash.load_neighborhood_minhash(infile)
    print('...done!')

    data = open(args.sigfile).read()
    siglist  = sourmash_signature.load_signatures(data, select_ksize=KSIZE)
    assert len(siglist) == 1
    sig = siglist[0]

    print(nbhd_mh.search(sig.estimator.mh, 0.2))
    

if __name__ == '__main__':
    main()
    #nbhd_mh, combined = main()
