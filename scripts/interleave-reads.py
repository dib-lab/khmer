#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Take two files containing left & right reads from a paired-end sequencing run,
and interleave them.

% scripts/interleave-reads.py <R1> <R2> [ -o <outputfile> ]

By default, output is sent to stdout; or use -o. Use '-h' for parameter help.
"""

# TODO: take fa as well?
#      support gzip option?

import screed
import sys
import itertools
import os
import argparse


def output_pair(r1, r2):
    if hasattr(r1, 'accuracy'):
        return "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n" % \
            (r1.name, r1.sequence, r1.accuracy,
             r2.name, r2.sequence, r2.accuracy)
    else:
        return ">%s\n%s\n>%s\n%s\n" % (r1.name, r1.sequence, r2.name,
                                       r2.sequence)


def main():
    parser = argparse.ArgumentParser(
        description='Produce interleaved files from R1/R2 paired files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('infiles', nargs='+')
    parser.add_argument('-o', '--output',
                        dest='output', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()

    s1_file = args.infiles[0]
    if len(args.infiles) == 2:
        s2_file = args.infiles[1]
    else:
        s2_file = s1_file.replace('_R1_', '_R2_')
        print >>sys.stderr, "given only one file;"
        " guessing that R2 file is %s" % s2_file

    fail = False
    if not os.path.exists(s1_file):
        print >>sys.stderr, "Error! R1 file %s does not exist" % s1_file
        fail = True

    if not os.path.exists(s2_file):
        print >>sys.stderr, "Error! R2 file %s does not exist" % s2_file
        fail = True

    if fail:
        sys.exit(1)

    print >>sys.stderr, "Interleaving:\n\t%s\n\t%s" % (s1_file, s2_file)

    n = 0
    for r1, r2 in itertools.izip(screed.open(s1_file), screed.open(s2_file)):
        if n % 100000 == 0:
            print >>sys.stderr, '...', n, 'pairs'
        n += 1

        name1 = r1.name
        if not name1.endswith('/1'):
            name1 += '/1'
        name2 = r2.name
        if not name2.endswith('/2'):
            name2 += '/2'

        assert name1[:-2] == name2[:-
                                   2], "This doesn't look like paired data!"
        " %s %s" % (name1, name2)

        r1.name = name1
        r2.name = name2
        args.output.write(output_pair(r1, r2))

    print >>sys.stderr, 'final: interleaved %d pairs' % n

if __name__ == '__main__':
    main()
