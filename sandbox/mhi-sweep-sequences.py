#! /usr/bin/env python
import sys
import os
import argparse
import khmer
import screed

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
    parser.add_argument('mhi_file')
    parser.add_argument('sigfile')
    parser.add_argument('seqfiles', nargs='+')
    args = parser.parse_args()

    print('loading nbhd minhashes...')
    nbhd_mh = khmer.load_neighborhood_minhash(args.mhi_file)
    print('...done!')
    
    data = open(args.sigfile).read()
    siglist = sourmash_signature.load_signatures(data, select_ksize=KSIZE)
    if not siglist:
        print('No signatures!')
        sys.exit(0)

    print('building ~chromosome level minhashes')
    combined = nbhd_mh.build_combined_minhashes(COMBINED_MH_SIZE)
    print('...built %d' % len(combined))

    # find appropriate MinHash (all above 0.9?)

    found = set()
    for c in combined:
        mh = c.get_minhash()
        for sig in siglist:
            if sig.estimator.mh.compare(mh) >= THRESHOLD:
                found.add(c)

    print('found %d combined MinHashes that match' % len(found))

    # extract tags

    tags = set()
    for c in found:
        t = c.get_tags()
        tags.update(t)

    print('a total of %d tags from %d combined MinHashes' % (len(tags),
                                                             len(found)))

    ng = khmer.Nodegraph(KSIZE, 1, 1)
    for t in tags:
        ng.add_tag(khmer.reverse_hash(t, KSIZE))

    # sweep through the seqfiles, extracting reads
    outfile = os.path.basename(args.sigfile) + '.mhi_sweep'
    outfp = open(outfile, 'w')
        
    total_extracted = 0
    total = 0
    for seqfile in args.seqfiles:
        m = 0
        for n, record in enumerate(screed.open(seqfile)):
            if n % 10000 == 0:
                print('...', n, m, seqfile, total_extracted)
                
            if ng.get_tags_and_positions(record.sequence):
                outfp.write('>%s\n%s\n' % (record.name, record.sequence))
                m += 1
        total_extracted += m
        total += n + 1

    print('Extracted %d sequences (out of %d total)'
              % (total_extracted, total))
    print('Extracted sequences are in %s' % outfile)
    
if __name__ == '__main__':
    main()
