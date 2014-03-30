#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Produce the k-mer abundance distribution for the given file.

% python scripts/abundance-dist.py [ -z -s ] <htname> <data> <histout>

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import khmer
import argparse
import os
from khmer.file import check_file_status, check_space


def main():
    parser = argparse.ArgumentParser(
        description="Output k-mer abundance distribution.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('hashname')
    parser.add_argument('datafile')
    parser.add_argument('histout')

    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')

    args = parser.parse_args()

    infiles = [args.hashname, args.datafile]
    for infile in infiles:
        check_file_status(infile)

    check_space(infiles)

    print('hashtable from', args.hashname)
    counting_hash = khmer.load_counting_hash(args.hashname)

    kmer_size = counting_hash.ksize()
    hashsizes = counting_hash.hashsizes()
    tracking = khmer._new_hashbits(kmer_size, hashsizes)

    print('K:', kmer_size)
    print('HT sizes:', hashsizes)
    print('outputting to', args.histout)

    if os.path.exists(args.histout):
        if not args.squash_output:
            print('ERROR: %s exists; not squashing.' % args.histout,
                  file=sys.stderr)
            sys.exit(1)

        print('** squashing existing file %s' % args.histout)

    print('preparing hist...')
    abundances = counting_hash.abundance_distribution(args.datafile, tracking)
    total = sum(abundances)

    if 0 == total:
        print("ERROR: abundance distribution is uniformly zero; "
              "nothing to report.", file=sys.stderr)
        print("\tPlease verify that the input files are valid.",
              file=sys.stderr)
        sys.exit(1)
    histout = args.histout
    hash_fp = open(histout, 'w')

    sofar = 0
    for _, i in enumerate(abundances):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print(_, i, sofar, round(frac, 3), file=hash_fp)

        if sofar == total:
            break

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
