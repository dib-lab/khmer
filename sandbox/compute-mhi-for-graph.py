#! /usr/bin/env python
import khmer
import screed
from cPickle import dump
import argparse
import os.path

import sys
sys.path.append('../sourmash')
import sourmash_lib, sourmash_signature

MH_LEN=50
KSIZE=32

def load_and_tag(ct, filename):
    print 'reading and tagging sequences'
    for record in screed.open(filename):
        print '.', record.name
        ct.consume_and_tag(record.sequence)
    print '...done!'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfile')
    args = parser.parse_args()
    
    seqfile = args.seqfile
    
    ct = khmer.Countgraph(KSIZE, 4e8, 2)
    ct._set_tag_density(200)

    ###

    load_and_tag(ct, seqfile)

    print('building nbhd & ~chromosome level minhashes')
    for n, mh in enumerate(ct.build_combined_minhashes()):
        e = sourmash_lib.Estimators(n=0, ksize=KSIZE)
        e.mh = mh
        sig = sourmash_signature.SourmashSignature('t@idyll.org', e)
        out = sourmash_signature.save_signatures([sig])
        open('%s.%d.sig' % (os.path.basename(seqfile), n + 1), 'w').write(out)
    print('wrote %d sigs' % (n + 1,))


if __name__ == '__main__':
    main()
