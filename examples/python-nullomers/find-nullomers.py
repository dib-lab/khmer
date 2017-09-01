#! /usr/bin/env python
"""
Find k-mers that have abundance 0 (fwd and rev comp) in input files.

Uses exact counting and requires 4**K / 8 bytes, so don't use this
for large k :).

Usage:

   ./find-nullomers.py -k 11 inp.fasta inp2.fasta inp3.fasta
"""
import sys
import khmer
import argparse
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('contigs', nargs='+')
    p.add_argument('-k', '--ksize', type=int, default=11)
    args = p.parse_args()

    assert args.ksize % 2 == 1, "K must be odd"
    assert args.ksize > 1, "K must be larger than 1"

    SIZE=4**args.ksize                               # use exact counting
    
    print('allocating lots of memory for exact counts: {} bytes'.format(SIZE / 8), file=sys.stderr)

    # use Nodegraph because we want to use the predictable hash function
    # that only goes up to 4**K.
    ct = khmer.Nodegraph(args.ksize, SIZE, 1)

    for filename in args.contigs:
        print("consuming '{}'".format(filename), file=sys.stderr)
        ct.consume_seqfile(filename)
    print('...done!', file=sys.stderr)

    print('Iterating over all {}-mers'.format(args.ksize), file=sys.stderr)

    # for large K, this is going to end up producing a massive amount of
    # output...
    for i in range(SIZE):
        # the point of this is to canonicalize k-mers to less of RC.
        kmer = ct.reverse_hash(i)
        if kmer.startswith('A') or kmer.startswith('C'):
            if ct.get(kmer) == 0:
                print(kmer)


if __name__ == '__main__':
    main()
