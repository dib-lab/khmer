#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org

from __future__ import print_function
import argparse
import screed
import sys
import khmer


def output_single(read):
    if hasattr(read, 'quality'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.quality)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-coverage', type=int, default=None)
    parser.add_argument('-M', '--max-coverage', type=int, default=None)
    parser.add_argument('input_counting_table')
    parser.add_argument('input_readfile')
    parser.add_argument('output_readfile')
    args = parser.parse_args()

    print('min_coverage: %s' % args.min_coverage, file=sys.stderr)
    print('max_coverage: %s' % args.max_coverage, file=sys.stderr)

    if not (args.min_coverage or args.max_coverage):
        print("neither min nor max coverage specified!? exiting!", file=sys.stderr)
        sys.exit(1)

    if args.min_coverage and args.max_coverage and \
       args.max_coverage < args.min_coverage:
        print("min_coverage > max_coverage!? exiting!", file=sys.stderr)
        sys.exit(1)

    htable = khmer.load_counting_hash(args.input_counting_table)
    output_file = args.output_readfile
    output_fp = open(output_file, 'w')

    n_kept = 0
    n = 0
    for n, record in enumerate(screed.open(args.input_readfile)):
        if n % 100000 == 0:
            print('...', n, n_kept, file=sys.stderr)

        seq = record.sequence.upper()
        seq = seq.replace('N', 'A')

        try:
            med, _, _ = htable.get_median_count(seq)
        except ValueError:
            continue

        keep = True
        if args.min_coverage and med < args.min_coverage:
            keep = False

        if args.max_coverage and med > args.max_coverage:
            keep = False

        if keep:
            n_kept += 1

            output_fp.write(output_single(record))

    print('consumed %d reads; kept %d' % (n, n_kept), file=sys.stderr)

if __name__ == '__main__':
    main()
