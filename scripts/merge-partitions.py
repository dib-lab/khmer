#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Merge multiple pmap files into a single one.

% python scripts/merge-partitions.py <base>

Load <base>.subset.*.pmap and merge into a single pmap file.  Final
merged pmap file will be in <base>.pmap.merged.
"""

import argparse
import glob
import os
import textwrap
import khmer
from khmer.file import check_file_status, check_space
from khmer.khmer_args import info

DEFAULT_K = 32


def get_parser():
    epilog = """
    Take the ${graphbase}.subset.#.pmap files and merge them all into a single
    ${graphbase}.pmap.merged file for :program:`annotate-partitions.py` to use.
    """
    parser = argparse.ArgumentParser(
        description="Merge partition map '.pmap' files.",
        epilog=textwrap.dedent(epilog))
    parser.add_argument('--ksize', '-k', type=int, default=DEFAULT_K,
                        help="k-mer size (default: %d)" % DEFAULT_K)
    parser.add_argument('--keep-subsets', dest='remove_subsets',
                        default=True, action='store_false',
                        help='Keep individual subsets (default: False)')
    parser.add_argument('graphbase', help='basename for input and output '
                        'files')
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    info('merge-partitions.py', ['graph'])
    args = get_parser().parse_args()

    output_file = args.graphbase + '.pmap.merged'
    pmap_files = glob.glob(args.graphbase + '.subset.*.pmap')

    print 'loading %d pmap files (first one: %s)' % (len(pmap_files),
                                                     pmap_files[0])

    ksize = args.ksize
    htable = khmer.new_hashbits(ksize, 1, 1)

    for _ in pmap_files:
        check_file_status(_)

    check_space(pmap_files)

    for pmap_file in pmap_files:
        print 'merging', pmap_file
        htable.merge_subset_from_disk(pmap_file)

    print 'saving merged to', output_file
    htable.save_partitionmap(output_file)

    if args.remove_subsets:
        print 'removing pmap files'
        for pmap_file in pmap_files:
            os.unlink(pmap_file)

if __name__ == '__main__':
    main()
