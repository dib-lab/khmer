#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Interleave left and right reads.

Take two files containing left & right reads from a paired-end sequencing run,
and interleave them.

% scripts/interleave-reads.py <R1> <R2> [ -o <outputfile> ]

By default, output is sent to stdout; or use -o. Use '-h' for parameter help.
"""
from __future__ import print_function

import screed
import sys
import os
import textwrap
import argparse
import khmer
from khmer.kfile import check_input_files, check_space, is_block
from khmer.khmer_args import info, sanitize_epilog
from khmer.kfile import (add_output_compression_type, get_file_writer,
                         describe_file_handle)
from khmer.utils import (write_record_pair, check_is_left, check_is_right,
                         check_is_pair)

try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest


def get_parser():
    epilog = """
    The output is an interleaved set of reads, with each read in <R1> paired
    with a read in <R2>. By default, the output goes to stdout unless
    :option:`-o`/:option:`--output` is specified.

    As a "bonus", this file ensures that if read names are not already
    formatted properly, they are reformatted consistently, such that
    they look like the pre-1.8 Casava format (@name/1, @name/2).

    Example::

             interleave-reads.py tests/test-data/paired.fq.1 \
                     tests/test-data/paired.fq.2 -o paired.fq"""
    parser = argparse.ArgumentParser(
        description='Produce interleaved files from R1/R2 paired files',
        epilog=textwrap.dedent(epilog),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('left')
    parser.add_argument('right')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=argparse.FileType('wb'),
                        default=sys.stdout)
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    add_output_compression_type(parser)
    return parser


def main():
    info('interleave-reads.py')
    args = sanitize_epilog(get_parser()).parse_args()

    check_input_files(args.left, args.force)
    check_input_files(args.right, args.force)
    check_space([args.left, args.right], args.force)

    s1_file = args.left
    s2_file = args.right

    fail = False
    if not os.path.exists(s1_file):
        print("Error! R1 file %s does not exist" % s1_file, file=sys.stderr)
        fail = True

    if not os.path.exists(s2_file):
        print("Error! R2 file %s does not exist" % s2_file, file=sys.stderr)
        fail = True

    if fail and not args.force:
        sys.exit(1)

    print("Interleaving:\n\t%s\n\t%s" % (s1_file, s2_file), file=sys.stderr)

    outfp = get_file_writer(args.output, args.gzip, args.bzip)

    counter = 0
    screed_iter_1 = screed.open(s1_file)
    screed_iter_2 = screed.open(s2_file)
    for read1, read2 in zip_longest(screed_iter_1, screed_iter_2):
        if read1 is None or read2 is None:
            print(("ERROR: Input files contain different number"
                   " of records."), file=sys.stderr)
            sys.exit(1)

        if counter % 100000 == 0:
            print('...', counter, 'pairs', file=sys.stderr)
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
            print("ERROR: This doesn't look like paired data! "
                  "%s %s" % (read1.name, read2.name), file=sys.stderr)
            sys.exit(1)

        write_record_pair(read1, read2, outfp)

    print('final: interleaved %d pairs' % counter, file=sys.stderr)
    print('output written to', describe_file_handle(outfp), file=sys.stderr)

if __name__ == '__main__':
    main()
