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
    parser.add_argument('-S', '--savegraph', action='store_true')
    parser.add_argument('--protein', action='store_true')
    parser.add_argument('-k', '--ksize', type=int, default=KSIZE)
    args = parser.parse_args()
    
    seqfile = args.seqfile
    outfile = args.output
    if not outfile:
        outfile = os.path.basename(args.seqfile) + '.mhi'

    print('loading sequences from', seqfile)
    print('will save MinHash index to', outfile)

    ct = khmer.Countgraph(args.ksize, 5e8, 2)
    ct._set_tag_density(200)

    ###

    load_and_tag(ct, seqfile)

    if args.savegraph:
        print('saving graph + tags...')
        ct.save(outfile + '.cg')
        ct.save_tagset(outfile + '.cg.tags')

    print('building nbhd minhashes...')
    nbhd_mh = ct.build_neighborhood_minhashes()   #@CTB args.protein
    nbhd_mh.save(outfile)
    print('...done building! mhi saved to', outfile)

if __name__ == '__main__':
    main()
