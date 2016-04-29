#! /usr/bin/env python
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
    raise
    pass

KSIZE=32
COMBINED_MH_SIZE=1000
TAGCOUNT_THRESHOLD=4000

def filter_combined(combined, min_tagcount=TAGCOUNT_THRESHOLD):
    cc = [ (len(c.get_tags()), c) for c in combined ]
    cc.sort(reverse=True, key=lambda x: x[0])
    print(cc[:10])
    return [ c for (n, c) in cc if n >= min_tagcount ]

def load_and_tag(ct, filename):
    print('reading and tagging sequences')
    for record in screed.open(filename):
        print('.', record.name)
        ct.consume_and_tag(record.sequence)
    print('...done!')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('-l', '--load-seqfile')
    args = parser.parse_args()
    
    infile = args.infile
    seqfile = args.load_seqfile

    print('will load MinHash index from', infile)

    print('loading nbhd minhashes...')
    nbhd_mh = khmer.load_neighborhood_minhash(infile)
    print('...done!')

    if seqfile:
        print('loading sequences from', seqfile)

        ct = khmer.Countgraph(KSIZE, 1, 1)
        ct._set_tag_density(200)

    ###

        load_and_tag(ct, seqfile)

        combined = []
        for record in screed.open(seqfile):
            print('.2', record.name)
            x = []
            for p, tag in ct.get_tags_and_positions(record.sequence):
                x.append(tag)

            combined.append(nbhd_mh.combine_from_tags(COMBINED_MH_SIZE, x))
        print(combined)
        basename = os.path.basename(seqfile)
    else:
        print('building ~chromosome level minhashes')
        combined = nbhd_mh.build_combined_minhashes(COMBINED_MH_SIZE)
        combined = filter_combined(combined)

        basename = os.path.basename(infile)
        if basename.endswith('.mhi'):
            basename = basename[:-4]
        

    for n, mh in enumerate(combined):
        e = sourmash_lib.Estimators(n=0, ksize=KSIZE)
        e.mh = mh.get_minhash()
        sig = sourmash_signature.SourmashSignature('t@idyll.org', e)
        out = sourmash_signature.save_signatures([sig])
        open('%s.%d.sig' % (basename, n + 1), 'w').write(out)
    if combined:
        print('wrote %d sigs to %s.%%d.sig' % (n + 1, basename))
    else:
        print("no signatures to extract")

if __name__ == '__main__':
    main()
