#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
"""
Calculate the mismatch error profile for shotgun data, using a subset of
reads.  The output is placed in <infile>.errhist in the cwd by default.

% sandbox/calc-error-profile.py [ -o outfile ] <infile>

Reads FASTQ and FASTA input.
"""

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
    ht = khmer.Countgraph(K, HASHSIZE, N_HT)

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
