#! /usr/bin/env python
from __future__ import print_function
import khmer
import screed
try:
    from cPickle import dump
except ImportError:
    from pickle import dump
import argparse
import os.path

import sys
sys.path.append('../sourmash')
try:
    import sourmash_lib, sourmash_signature
except ImportError:
    pass

KSIZE=32
COMBINED_MH_SIZE=500
COMBINE_THIS_MANY=10000

def load_and_tag(ct, filename):
    print('reading and tagging sequences')
    for record in screed.open(filename):
        print('.', record.name)
        ct.consume_and_tag(record.sequence)
    print('...done!')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('sigfile')
    parser.add_argument('-S', '--combine-this-many', type=int,
                        default=COMBINE_THIS_MANY)
    args = parser.parse_args()
    
    infile = args.infile
    
    sigdata = open(args.sigfile).read()
    sigs = sourmash_signature.load_signatures(sigdata, select_ksize=KSIZE)
    assert len(sigs) == 1
    sig = sigs[0]

    print('will load MinHash index from', infile)

    print('loading nbhd minhashes...')
    nbhd_mh = khmer.load_neighborhood_minhash(infile)
    print('...done!')

    print('building ~chromosome level minhashes')
    combined = nbhd_mh.build_combined_minhashes(args.combine_this_many,
                                                COMBINED_MH_SIZE)

    filtered = []
    mh = sig.estimator.mh
    for x in combined:
        cmh = x.get_minhash()
        if mh.compare(cmh) > 0.01:
            filtered.append(x)

    print('found %d of %d' % (len(filtered), len(combined)))

    tagset = set()
    for x in filtered:
        tagset.update(x.get_tags())

    merged = nbhd_mh.combine_from_tags(2*COMBINED_MH_SIZE, list(tagset))

    basename = os.path.basename(infile)
    if basename.endswith('.mhi'):
        basename = basename[:-4]

    print(merged.get_minhash().compare(sig.estimator.mh))
    print(sig.estimator.mh.compare(merged.get_minhash()))

    return nbhd_mh, combined

if __name__ == '__main__':
    nbhd_mh, combined = main()
