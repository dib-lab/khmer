#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Annotate sequences with partition numbers.

% python scripts/annotate-partitions.py <pmap_file> <file1> [ <file2> ... ]

Partition-annotated sequences will be in <fileN>.part.

Use '-h' for parameter help.
"""
from __future__ import print_function

import os
import argparse
import textwrap
import khmer
import sys
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import info

DEFAULT_K = 32


def get_parser():
    epilog = """
    Load in a partitionmap (generally produced by partition-graph.py or
    merge-partitions.py) and annotate the sequences in the given files with
    their partition IDs. Use :program:`extract-partitions.py` to extract
    sequences into separate group files.

    Example (results will be in ``random-20-a.fa.part``)::

        load-graph.py -k 20 example tests/test-data/random-20-a.fa
        partition-graph.py example
        merge-partitions.py -k 20 example
        annotate-partitions.py -k 20 example tests/test-data/random-20-a.fa
    """
    parser = argparse.ArgumentParser(
        description="Annotate sequences with partition IDs.",
        epilog=textwrap.dedent(epilog),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--ksize', '-k', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('graphbase', help='basename for input and output '
                        'files')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+', help='input FAST[AQ] sequences to '
                        'annotate.')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                        khmer.__version__)
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('annotate-partitions.py', ['graph'])
    args = get_parser().parse_args()

    ksize = args.ksize
    filenames = args.input_filenames
    htable = khmer.new_hashbits(ksize, 1, 1)

    partitionmap_file = args.graphbase + '.pmap.merged'

    check_input_files(partitionmap_file, args.force)
    for _ in filenames:
        check_input_files(_, args.force)

    check_space(filenames, args.force)

    print('loading partition map from:', partitionmap_file, file=sys.stderr)
    htable.load_partitionmap(partitionmap_file)

    for infile in filenames:
        print('outputting partitions for', infile, file=sys.stderr)
        outfile = os.path.basename(infile) + '.part'
        part_count = htable.output_partitions(infile, outfile)
        print('output %d partitions for %s' % (
            part_count, infile), file=sys.stderr)
        print('partitions are in', outfile, file=sys.stderr)

if __name__ == '__main__':
    main()
