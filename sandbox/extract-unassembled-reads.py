#! /usr/bin/env python
"""
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

    # first, consume & tag the reads

    for readfile in args.readfiles:
        print('loading & tagging reads from:', readfile)
        ng.consume_seqfile_and_tag(readfile)

    ## next, consume & label the assembly

    print('loading & tagging assembly from:', args.assembly)
    lh = khmer._GraphLabels(ng)
    lh.consume_seqfile_and_tag_with_labels(args.assembly)

    if args.output:
        outfp = open(args.output, 'w')

    ## finally, walk across the reads & find those with no labels
        
    n = 0
    m = 0

    for readfile in args.readfiles:
        print('loading reads from:', readfile)
        if not args.output:
            outfile = os.path.basename(readfile) + '.leftover2'
            outfp = open(outfile, 'w')
            print('writing to:', outfile, file=sys.stderr)

        for record in screed.open(readfile):
            if n % 100000 == 0 and n:
                print('...', readfile, n, m, file=sys.stderr)
            x = ng.get_tags_and_positions(record.sequence)

            do_extract = False
            for (pos, tag) in x:
                if not lh.get_tag_labels(tag):
                    do_extract = True
                    break

            if do_extract:
                khmer.utils.write_record(record, outfp)
                m += 1

            n += 1

        if not args.output:
            outfp.close()

    print('%d left out of assembly, of %d reads' % (m, n), file=sys.stderr)


if __name__ == '__main__':
    main()
