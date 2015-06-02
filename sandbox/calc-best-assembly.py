#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE. Contact: ctb@msu.edu
#
from __future__ import print_function
import screed
import argparse
import sys

DEFAULT_SIZE_CUTOFF=500

def calculate_bp_above_cutoff(filename, cutoff):
    total = 0
    for record in screed.open(filename):
        if len(record.sequence) >= cutoff:
            total += len(record.sequence)
    return total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_SIZE_CUTOFF)
    parser.add_argument('-o', '--output-file', dest='output_file',
                        type=argparse.FileType('w'))
    parser.add_argument('-q', '--quiet', dest='quiet',
                        type=bool)
    parser.add_argument('assembly_files', nargs='+')

    args = parser.parse_args()

    stats = []
    for filename in args.assembly_files:
        try:
            total = calculate_bp_above_cutoff(filename, args.cutoff)
        except IOError:
            print("** WARNING: %s does not exist, skipping" %\
                filename, file=sys.stderr)
            continue

        stats.append((total, filename))

        if not args.quiet:
            print("assembly %s has %d bp > %d" % (filename,
                                                                total,
                                                                args.cutoff), file=sys.stderr)

    stats.sort(reverse=True)

    best_total, winner_file = stats[0]
    print('----', file=sys.stderr)
    print("assembly %s wins: %d total bp > %d" % (winner_file,
                                                                best_total,
                                                                args.cutoff), file=sys.stderr)

    if args.output_file:
        for record in screed.open(winner_file, parse_description=False):
            print('>%s\n%s' % (record.name,
                                                   record.sequence), file=args.output_file)

    print(winner_file)

if __name__ == '__main__':
    main()
