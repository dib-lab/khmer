#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
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
            print >>sys.stderr, "** WARNING: %s does not exist, skipping" %\
                filename
            continue

        stats.append((total, filename))

        if not args.quiet:
            print >>sys.stderr, "assembly %s has %d bp > %d" % (filename,
                                                                total,
                                                                args.cutoff)

    stats.sort(reverse=True)

    best_total, winner_file = stats[0]
    print >>sys.stderr, '----'
    print >>sys.stderr, "assembly %s wins: %d total bp > %d" % (winner_file,
                                                                best_total,
                                                                args.cutoff)

    if args.output_file:
        for record in screed.open(winner_file, parse_description=False):
            print >>args.output_file, '>%s\n%s' % (record.name,
                                                   record.sequence)

    print winner_file

if __name__ == '__main__':
    main()
