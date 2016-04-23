#! /usr/bin/env python
from __future__ import print_function
import khmer
import screed
from cPickle import dump
import argparse
import os.path

import sys
sys.path.append('../sourmash')
try:
    import sourmash_lib, sourmash_signature
except ImportError:
    pass

MH_LEN=50
KSIZE=32

def load_and_tag(ct, filename):
    print('reading and tagging sequences')
    for record in screed.open(filename):
        print('.', record.name)
        ct.consume_and_tag(record.sequence)
    print('...done loading sequences!')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfile')
    parser.add_argument('-o', '--output', metavar='output_filename')
    args = parser.parse_args()
    
    seqfile = args.seqfile
    outfile = args.output
    if not outfile:
        outfile = os.path.basename(args.seqfile) + '.mhi'

    print('loading sequences from', seqfile)
    print('will save MinHash index to', outfile)

    ct = khmer.Countgraph(KSIZE, 4e8, 2)
    ct._set_tag_density(200)

    ###

    load_and_tag(ct, seqfile)

    print('building nbhd minhashes...')
    nbhd_mh = ct.build_neighborhood_minhashes(20, 9999999967)
    nbhd_mh.save(outfile)
    print('...done building! mhi saved to', outfile)

if __name__ == '__main__':
    main()
