#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
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
from khmer.file import check_file_status, check_space


def output_pair(read1, read2):
    if hasattr(read1, 'accuracy'):
        return "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n" % \
            (read1.name, read1.sequence, read1.accuracy,
             read2.name, read2.sequence, read2.accuracy)
    else:
        return ">%s\n%s\n>%s\n%s\n" % (read1.name, read1.sequence, read2.name,
                                       read2.sequence)


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
    return parser


def main():
    args = get_parser().parse_args()

    for _ in args.infiles:
        check_file_status(_)

    check_space(args.infiles)

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

    if fail:
        sys.exit(1)

    print >> sys.stderr, "Interleaving:\n\t%s\n\t%s" % (s1_file, s2_file)

    counter = 0
    for read1, read2 in itertools.izip(screed.open(s1_file),
                                       screed.open(s2_file)):
        if counter % 100000 == 0:
            print >> sys.stderr, '...', counter, 'pairs'
        counter += 1

        name1 = read1.name
        if not name1.endswith('/1'):
            name1 += '/1'
        name2 = read2.name
        if not name2.endswith('/2'):
            name2 += '/2'

        assert name1[:-2] == name2[:-2], \
            "This doesn't look like paired data! %s %s" % (name1, name2)

        read1.name = name1
        read2.name = name2
        args.output.write(output_pair(read1, read2))

    print >> sys.stderr, 'final: interleaved %d pairs' % counter

if __name__ == '__main__':
    main()
