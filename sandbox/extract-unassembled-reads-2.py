#! /usr/bin/env python
"""
This script differs from extract-unassembled-reads.py in that it only tags
the assembly, not the reads.  This is more efficient (and streaming for
the reads!) but it will erroneously extract some small fraction of reads
because they miss the tags for reasons of length or errors.

Procedure:

* hard trim the reads at an abundance of ~5 (to avoid low abundance, and
    eliminate erroneous paths)
* and/or variable-coverage trim the reads at an abundance of 3, to eliminate
    erroneous paths from super-high-abundance data
* run this script with the assembly & the remaining reads.
"""
import sys
import os.path
import khmer, khmer.utils
import screed
import argparse


DEFAULT_KSIZE=31
NODEGRAPH_SIZE=1e8


def main():
    p = argparse.ArgumentParser()
    p.add_argument('assembly')
    p.add_argument('readfiles', nargs='+')
    p.add_argument('-o', '--output', default=None)
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('-x', '--tablesize', default=NODEGRAPH_SIZE,
                   type=float)
    args = p.parse_args()

    ng = khmer.Nodegraph(args.ksize, args.tablesize, 4)
    ng._set_tag_density(20)

    print('loading & tagging assembly from:', args.assembly)
    ng.consume_seqfile_and_tag(args.assembly)

    if args.output:
        outfp = open(args.output, 'w')

    n = 0
    m = 0
    for readfile in args.readfiles:
        print('loading reads from:', readfile)
        if not args.output:
            outfile = os.path.basename(readfile) + '.leftover'
            outfp = open(outfile, 'w')
            print('writing to:', outfile, file=sys.stderr)

        for record in screed.open(readfile):
            if n % 100000 == 0 and n:
                print('...', readfile, n, m, file=sys.stderr)
            x = ng.get_tags_and_positions(record.sequence)
            if not x:
                khmer.utils.write_record(record, outfp)
                m += 1
            n += 1

        if not args.output:
            outfp.close()

    print('%d left out of assembly, of %d reads' % (m, n), file=sys.stderr)


if __name__ == '__main__':
    main()
