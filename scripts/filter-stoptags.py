#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Trim sequences at k-mers in the given stoptags file.  Output sequences
will be placed in 'infile.stopfilt'.

% python scripts/filter-stoptags.py <stoptags> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import os
import khmer
import argparse
import textwrap
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer.file import check_file_status, check_space
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
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    info('filter-stoptags.py', ['graph'])
    args = get_parser().parse_args()
    stoptags = args.stoptags_file
    infiles = args.input_filenames

    for _ in infiles:
        check_file_status(_)

    check_space(infiles)

    print 'loading stop tags, with K', args.ksize
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
        print 'filtering', infile
        outfile = os.path.basename(infile) + '.stopfilt'

        outfp = open(outfile, 'w')

        tsp = ThreadedSequenceProcessor(process_fn)
        tsp.start(verbose_loader(infile), outfp)

        print 'output in', outfile

if __name__ == '__main__':
    main()
