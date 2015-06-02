#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Sequence trimming using stoptags.

Trim sequences at k-mers in the given stoptags file.  Output sequences
will be placed in 'infile.stopfilt'.

% python scripts/filter-stoptags.py <stoptags> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
from __future__ import print_function

import os
import khmer
import argparse
import textwrap
import sys
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info

# @CTB K should be loaded from file...
DEFAULT_K = 32


def get_parser():
    epilog = """
    Load stoptags in from the given .stoptags file and use them to trim
    or remove the sequences in <file1-N>.  Trimmed sequences will be placed in
    <fileN>.stopfilt.
    """
    parser = argparse.ArgumentParser(
        description="Trim sequences at stoptags.",
        epilog=textwrap.dedent(epilog),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ksize', '-k', default=DEFAULT_K, type=int,
                        help='k-mer size')
    parser.add_argument('stoptags_file', metavar='input_stoptags_filename')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('filter-stoptags.py', ['graph'])
    args = get_parser().parse_args()
    stoptags = args.stoptags_file
    infiles = args.input_filenames

    for _ in infiles:
        check_input_files(_, args.force)

    check_space(infiles, args.force)

    print('loading stop tags, with K', args.ksize, file=sys.stderr)
    htable = khmer.new_hashbits(args.ksize, 1, 1)
    htable.load_stop_tags(stoptags)

    def process_fn(record):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = htable.trim_on_stoptags(seq)

        if trim_at >= args.ksize:
            return name, trim_seq

        return None, None

    # the filtering loop
    for infile in infiles:
        print('filtering', infile, file=sys.stderr)
        outfile = os.path.basename(infile) + '.stopfilt'

        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), outfp)

        print('output in', outfile, file=sys.stderr)

if __name__ == '__main__':
    main()
