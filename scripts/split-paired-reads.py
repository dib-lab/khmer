#! /usr/bin/env python
#
# This script is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
De-interleave a file.

Take an interleaved set of reads (/1 and /2), and extract them into separate
files (.1 and .2).

% scripts/split-paired-reads.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""
from __future__ import print_function
import screed
import sys
import os
import textwrap
import argparse
import khmer
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info
from khmer.utils import (write_record, check_is_left, check_is_right,
                         broken_paired_reader)


def get_parser():
    epilog = """
    Some programs want paired-end read input in the One True Format, which is
    interleaved; other programs want input in the Insanely Bad Format, with
    left- and right- reads separated. This reformats the former to the latter.

    The directory into which the left- and right- reads are output may be
    specified using :option:`-o`/:option:`--output-dir`. This directory will be
    created if it does not already exist.

    Alternatively, you can specify the filenames directly with
    :option:`-1`/:option:`--output-first` and
    :option:`-2`/:option:`--output-second`, which will override the
    :option:`-o`/:option:`--output-dir` setting on a file-specific basis.

    :option:`-p`/:option:`--force-paired` will require the input file to
    be properly interleaved; by default, this is not required.

    :option:`-0`/:option:'`--output-orphans` will separate orphans into 
    the specified file

    Example::

        split-paired-reads.py tests/test-data/paired.fq

    Example::

        split-paired-reads.py -o ~/reads-go-here tests/test-data/paired.fq

    Example::

        split-paired-reads.py -1 reads.1 -2 reads.2 tests/test-data/paired.fq
    """
    parser = argparse.ArgumentParser(
        description='Split interleaved reads into two files, left and right.',
        epilog=textwrap.dedent(epilog),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('infile', nargs='?', default='/dev/stdin')

    parser.add_argument('-o', '--output-dir', metavar="output_directory",
                        dest='output_directory', default='', help='Output '
                        'split reads to specified directory. Creates '
                        'directory if necessary')
    parser.add_argument('-0', '--output-orphaned', metavar='output_orphaned',
                        nargs='?', default=False, const=True, 
                        help='Allow "orphaned" reads and extract them to ' + \
                        'this file',
                        type=argparse.FileType('wb'))
    parser.add_argument('-1', '--output-first', metavar='output_first',
                        default=None, help='Output "left" reads to this '
                        'file', type=argparse.FileType('wb'))
    parser.add_argument('-2', '--output-second', metavar='output_second',
                        default=None, help='Output "right" reads to this '
                        'file', type=argparse.FileType('wb'))

    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('split-paired-reads.py')
    args = get_parser().parse_args()

    infile = args.infile

    filenames = [infile]
    check_input_files(infile, args.force)
    check_space(filenames, args.force)
    out0 = None
    fp_out0 = None

    basename = os.path.basename(infile)

    # decide where to put output files - specific directory? or just default?
    if infile in ('/dev/stdin', '-'):
        if not (args.output_first and args.output_second) or \
           (args.output_orphaned is not True
            and args.output_orphaned is not False):
            print("Accepting input from stdin; "
                  "output filenames must be provided.", file=sys.stderr)
            sys.exit(1)
    elif args.output_directory:
        if not os.path.exists(args.output_directory):
            os.makedirs(args.output_directory)
        out1 = os.path.join(args.output_directory, basename + '.1')
        out2 = os.path.join(args.output_directory, basename + '.2')
        if args.output_orphaned is True:
            out0 = os.path.join(args.output_directory, basename + '.0')
    else:
        out1 = basename + '.1'
        out2 = basename + '.2'
        if args.output_orphaned is True:
            out0 = basename + '.0'

    # OVERRIDE output file locations with -1, -2
    if args.output_first:
        fp_out1 = args.output_first
        out1 = fp_out1.name
    else:
        # Use default filename created above
        fp_out1 = open(out1, 'w')
    if args.output_second:
        fp_out2 = args.output_second
        out2 = fp_out2.name
    else:
        # Use default filename created above
        fp_out2 = open(out2, 'w')

    # put orphaned reads here!
    allow_orphans = False
    if args.output_orphaned is not True and args.output_orphaned is not False:
        fp_out0 = args.output_orphaned
        out0 = fp_out0.name
        allow_orphans = True
    elif args.output_orphaned is True:
        fp_out0 = open(out0, 'w')
        allow_orphans = True

    counter1 = 0
    counter2 = 0
    counter3 = 0
    index = None

    screed_iter = screed.open(infile)

    # walk through all the reads in broken-paired mode.
    paired_iter = broken_paired_reader(screed_iter)
    for index, is_pair, record1, record2 in paired_iter:
        if index % 10000 == 0:
            print('...', index, file=sys.stderr)

        if is_pair:
            write_record(record1, fp_out1)
            counter1 += 1
            write_record(record2, fp_out2)
            counter2 += 1
        elif allow_orphans:
            write_record(record1, fp_out0)
            counter3 += 1
        else:
            print('ERROR, %s is not part of a pair' %
                  record1.name, file=sys.stderr)
            sys.exit(1)

    print("DONE; split %d sequences (%d left, %d right, %d orphans)" %
          (counter1 + counter2, counter1, counter2, counter3), file=sys.stderr)
    print("/1 reads in %s" % out1, file=sys.stderr)
    print("/2 reads in %s" % out2, file=sys.stderr)
    if allow_orphans:
        print("orphans in %s" % out0, file=sys.stderr)

if __name__ == '__main__':
    main()
