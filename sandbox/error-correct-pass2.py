#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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
Error correct reads based on a counting hash from a diginorm step.
Output sequences will be put in inputfile.corr.

% python scripts/error-correct-pass2 <counting.ct> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
import sys
import os
import screed
import khmer
from khmer import Countgraph
from khmer import khmer_args
from khmer.khmer_args import FileType as khFileType

DEFAULT_CUTOFF = 2


def output_single(read, new_sequence):
    name = read.name
    sequence = new_sequence

    quality = None
    if hasattr(read, 'quality'):
        quality = read.quality[:len(sequence)]
        sequence = sequence[:len(quality)]  # sequence is _lengthened_

    if quality:
        assert len(sequence) == len(quality), (sequence, quality)
        return "@%s\n%s\n+\n%s\n" % (name, sequence, quality)
    else:
        return ">%s\n%s\n" % (name, sequence)


def main():
    parser = khmer_args.build_counting_args(
        "Correct reads against an already-computed table",
        citations=['counting', 'SeqAn'])

    parser.add_argument("--trusted-cov", dest="trusted_cov", type=int,
                        default=DEFAULT_CUTOFF)
    parser.add_argument("--theta", dest="bits_theta", type=float, default=1.0)
    parser.add_argument('-o', '--output', dest='output_file',
                        help="output file for histogram; defaults to "
                             "<first filename>.corr in cwd.",
                        type=khFileType('w'), default=None)

    parser.add_argument('counts_table')
    parser.add_argument('readfile')

    args = parser.parse_args()

    print('loading counts')
    ht = Countgraph.load(args.counts_table)

    aligner = khmer.ReadAligner(ht,
                                args.trusted_cov,
                                args.bits_theta)

    print("trusted:", args.trusted_cov)

    corrfp = args.output_file
    if not corrfp:
        outfile = os.path.basename(args.readfile) + '.corr'
        corrfp = open(outfile, 'w')

    n_corrected = 0
    for n, read in enumerate(screed.open(args.readfile)):
        if n % 10000 == 0:
            print('...', n, n_corrected, file=sys.stderr)
        seq = read.sequence.replace('N', 'A')

        # build the alignment...
        score, graph_alignment, read_alignment, truncated = \
            aligner.align(seq)

        if not truncated:
            graph_seq = graph_alignment.replace("-", "")
            if graph_seq != seq:
                n_corrected += 1

            seq = graph_seq

        corrfp.write(output_single(read, seq))


if __name__ == '__main__':
    main()
