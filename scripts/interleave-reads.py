#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
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
import textwrap
import argparse
import khmer
from khmer.kfile import check_file_status, check_space
from khmer.khmer_args import info
from khmer.utils import (write_record_pair, check_is_left, check_is_right,
                         check_is_pair)


def get_parser():
    epilog = """
    The output is an interleaved set of reads, with each read in <R1> paired
    with a read in <R2>. By default, the output goes to stdout unless
    :option:`-o`/:option:`--output` is specified.

    As a "bonus", this file ensures that read names are formatted in a
    consistent way, such that they look like the pre-1.8 Casava format
    (@name/1, @name/2).

    Example::

""" "        interleave-reads.py tests/test-data/paired.fq.1 tests/test-data/paired.fq.2 -o paired.fq"  # noqa
    parser = argparse.ArgumentParser(
        description='Produce interleaved files from R1/R2 paired files',
        epilog=textwrap.dedent(epilog),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('infiles', nargs='+')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('interleave-reads.py')
    args = get_parser().parse_args()

    for _ in args.infiles:
        check_file_status(_, args.force)

    check_space(args.infiles, args.force)

    s1_file = args.infiles[0]
    if len(args.infiles) == 2:
        s2_file = args.infiles[1]
    else:
        s2_file = s1_file.replace('_R1_', '_R2_')
        print >> sys.stderr, ("given only one file; "
                              "guessing that R2 file is %s" % s2_file)

    fail = False
    if not os.path.exists(s1_file):
        print >> sys.stderr, "Error! R1 file %s does not exist" % s1_file
        fail = True

    if not os.path.exists(s2_file):
        print >> sys.stderr, "Error! R2 file %s does not exist" % s2_file
        fail = True

    if fail and not args.force:
        sys.exit(1)

    print >> sys.stderr, "Interleaving:\n\t%s\n\t%s" % (s1_file, s2_file)

    counter = 0
    screed_iter_1 = screed.open(s1_file, parse_description=False)
    screed_iter_2 = screed.open(s2_file, parse_description=False)
    for read1, read2 in itertools.izip(screed_iter_1, screed_iter_2):
        if counter % 100000 == 0:
            print >> sys.stderr, '...', counter, 'pairs'
        counter += 1

        name1 = read1.name
        if not check_is_left(name1):
            name1 += '/1'
        name2 = read2.name
        if not check_is_right(name2):
            name2 += '/2'

        read1.name = name1
        read2.name = name2

        if not check_is_pair(read1, read2):
            print >>sys.stderr, "This doesn't look like paired data! %s %s" % \
                (read1.name, read2.name)
            sys.exit(1)

        write_record_pair(read1, read2, args.output)

    print >> sys.stderr, 'final: interleaved %d pairs' % counter
    print >> sys.stderr, 'output written to', args.output

if __name__ == '__main__':
    main()
