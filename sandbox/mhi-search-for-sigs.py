#! /usr/bin/env python
"search an mhi for a whole bunch of .sigs"
from __future__ import print_function
import khmer
import screed
import argparse
import os.path

import sys
sys.path.append('../sourmash')
try:
    import sourmash_lib, sourmash_signature
except ImportError:
    pass

KSIZE=32
COMBINED_MH_SIZE=1000
THRESHOLD=0.1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('minhash_index')
    parser.add_argument('sigfiles', nargs='+')
    args = parser.parse_args()
    
    infile = args.minhash_index
    
    print('loading nbhd minhashes...')
    nbhd_mh = khmer._minhash.load_neighborhood_minhash(infile)
    print('...done!')
    
    print('building ~chromosome level minhashes')
    combined = nbhd_mh.build_combined_minhashes2(COMBINED_MH_SIZE)

    siglist = []
    for sigfile in args.sigfiles:
        data = open(sigfile).read()
        siglist.extend(sourmash_signature.load_signatures(data,
                                                          select_ksize=KSIZE))

    results = set()
    for sig in siglist:
        mh = sig.estimator.mh
        for x in combined:
            graph_mh = x.get_minhash()
            if graph_mh.compare(mh) > THRESHOLD:
                results.add(sig)

    print('found match to following signatures at threshold', THRESHOLD)
    for r in results:
        print('\t', r.name())


if __name__ == '__main__':
    main()
    #nbhd_mh, combined = main()
