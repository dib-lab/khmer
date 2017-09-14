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
import screed
import argparse
import sys

DEFAULT_SIZE_CUTOFF=500

def calculate_bp_above_cutoff(filename, cutoff):
    total = 0
    for record in screed.open(filename):
        if len(record.sequence) >= cutoff:
            total += len(record.sequence)
    return total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-C', '--cutoff', type=int, dest='cutoff',
                        default=DEFAULT_SIZE_CUTOFF)
    parser.add_argument('-o', '--output-file', dest='output_file',
                        type=argparse.FileType('w'))
    parser.add_argument('-q', '--quiet', dest='quiet',
                        type=bool)
    parser.add_argument('assembly_files', nargs='+')

    args = parser.parse_args()

    stats = []
    for filename in args.assembly_files:
        try:
            total = calculate_bp_above_cutoff(filename, args.cutoff)
        except IOError:
            print("** WARNING: %s does not exist, skipping" %\
                filename, file=sys.stderr)
            continue

        stats.append((total, filename))

        if not args.quiet:
            print("assembly %s has %d bp > %d" % (filename, total,
                                                  args.cutoff),
                  file=sys.stderr)

    stats.sort(reverse=True)

    best_total, winner_file = stats[0]
    print('----', file=sys.stderr)
    print("assembly %s wins: %d total bp > %d" % (winner_file,
                                                  best_total,
                                                  args.cutoff),
          file=sys.stderr)

    if args.output_file:
        for record in screed.open(winner_file):
            print('>%s\n%s' % (record.name, record.sequence),
                  file=args.output_file)

    print(winner_file)

if __name__ == '__main__':
    main()
