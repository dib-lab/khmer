#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014-2015, Michigan State University.
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
import argparse
import screed
import sys
import khmer


def output_single(read):
    if hasattr(read, 'quality'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.quality)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-coverage', type=int, default=None)
    parser.add_argument('-M', '--max-coverage', type=int, default=None)
    parser.add_argument('input_count_graph')
    parser.add_argument('input_readfile')
    parser.add_argument('output_readfile')
    args = parser.parse_args()

    print('min_coverage: %s' % args.min_coverage, file=sys.stderr)
    print('max_coverage: %s' % args.max_coverage, file=sys.stderr)

    if not (args.min_coverage or args.max_coverage):
        print("neither min nor max coverage specified!? exiting!", file=sys.stderr)
        sys.exit(1)

    if args.min_coverage and args.max_coverage and \
       args.max_coverage < args.min_coverage:
        print("min_coverage > max_coverage!? exiting!", file=sys.stderr)
        sys.exit(1)

    htable = khmer.load_countgraph(args.input_count_graph)
    output_file = args.output_readfile
    output_fp = open(output_file, 'w')

    n_kept = 0
    n = 0
    for n, record in enumerate(screed.open(args.input_readfile)):
        if n % 100000 == 0:
            print('...', n, n_kept, file=sys.stderr)

        seq = record.sequence.upper()
        seq = seq.replace('N', 'A')

        try:
            med, _, _ = htable.get_median_count(seq)
        except ValueError:
            continue

        keep = True
        if args.min_coverage and med < args.min_coverage:
            keep = False

        if args.max_coverage and med > args.max_coverage:
            keep = False

        if keep:
            n_kept += 1

            output_fp.write(output_single(record))

    print('consumed %d reads; kept %d' % (n, n_kept), file=sys.stderr)

if __name__ == '__main__':
    main()
