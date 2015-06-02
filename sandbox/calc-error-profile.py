#! /usr/bin/env python
#
# This script is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE. Contact: ctb@msu.edu
#
"""
Calculate the mismatch error profile for shotgun data, using a subset of
reads.  The output is placed in <infile>.errhist in the cwd by default.

% scripts/calc-error-profile.py [ -o outfile ] <infile>

Reads FASTQ and FASTA input.
"""
from __future__ import division
from __future__ import print_function

import sys
import argparse
import khmer
import screed
import os.path

N_HT = 4
HASHSIZE = 1e7
K = 20
C = 10

MAX_SEQ_LEN = 65535
MAX_READS = 1e8
CHECK_EXIT = 25000


def exit_condition(n_consumed, n_checked):
    return (n_checked >= n_consumed or
            n_checked > 2e5)


def main():
    parser = argparse.ArgumentParser(
        "Calculate read error profile based on k-mer "
        "abundances of shotgun data.")

    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', dest='output_file',
                        help="output file for histogram; defaults to "
                             "<first filename>.errhist in cwd.",
                        type=argparse.FileType('w'), default=None)
    parser.add_argument('--errors-per-read', dest='errors_per_read',
                        type=argparse.FileType('w'), default=None)

    args = parser.parse_args()

    #
    # Figure out what the output filename is going to be
    #

    output_file = args.output_file
    if output_file:
        output_filename = output_file.name
    else:
        filename = args.filenames[0]
        output_filename = os.path.basename(filename) + '.errhist'
        output_file = open(output_filename, 'w')

    # Start!

    # build a small counting hash w/default parameters. In general there
    # should be no need to change these parameters.
    ht = khmer.new_counting_hash(K, HASHSIZE, N_HT)

    # initialize list to contain counts of errors by position
    positions = [0] * MAX_SEQ_LEN
    lengths = []                  # keep track of sequence lengths

    n_consumed = 0
    bp_consumed = 0
    total = 0
    n_checked = 0

    # run through all the files; pick out reads; once they saturate,
    # look for errors.
    total = 0
    for filename in args.filenames:
        print('opening', filename, file=sys.stderr)

        for n, record in enumerate(screed.open(filename)):
            total += 1

            if total % CHECK_EXIT == 0:
                print('...', total, n_consumed, n_checked, file=sys.stderr)

                # two exit conditions: first, have we hit our max reads limit?
                if total >= MAX_READS:
                    break

                # OR, alternatively, have we counted enough reads?
                if exit_condition(n_consumed, n_checked):
                    break

            # for each sequence, calculate its coverage:
            seq = record.sequence.replace('N', 'A')
            med, _, _ = ht.get_median_count(seq)

            # if the coverage is unsaturated, consume.
            if med < C:
                ht.consume(seq)
                n_consumed += 1
                bp_consumed += len(seq)
            else:
                # for saturated data, find low-abund k-mers
                posns = ht.find_spectral_error_positions(seq, 2)
                lengths.append(len(seq))

                if args.errors_per_read:
                    print(record.name, \
                        ",".join(map(str, posns)), file=args.errors_per_read)

                # track the positions => errors
                for p in posns:
                    positions[p] += 1

                n_checked += 1

    # normalize for length
    lengths.sort()
    max_length = lengths[-1]

    length_count = [0] * max_length
    for j in range(max_length):
        length_count[j] = sum([1 for i in lengths if i >= j])

    # write!
    output_file.write('position error_count error_fraction\n')
    for n, i in enumerate(positions[:max_length]):
        print(n, i, float(i) / float(length_count[n]), file=output_file)

    output_file.close()

    print('', file=sys.stderr)
    print('total sequences:', total, file=sys.stderr)
    print('n consumed:', n_consumed, file=sys.stderr)
    print('n checked:', n_checked, file=sys.stderr)
    print('bp consumed:', bp_consumed, bp_consumed / float(C), file=sys.stderr)
    print('error rate: %.2f%%' % \
        (100.0 * sum(positions) / float(sum(lengths))), file=sys.stderr)

    print('Error histogram is in %s' % output_filename, file=sys.stderr)

    if not exit_condition(n_consumed, n_checked):
        print("", file=sys.stderr)
        print("** WARNING: not enough reads to get a good result", file=sys.stderr)
        print("** Is this high diversity sample / small subset?", file=sys.stderr)
        sys.exit(-1)


if __name__ == '__main__':
    main()
