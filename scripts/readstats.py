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
# pylint: disable=missing-docstring
"""
Display summary statistics for one or more FASTA/FASTQ files.

% scripts/readstats.py [ -o output.txt ] <file1> <file2>

Use '-h' for parameter help.
"""

import argparse
import sys
import csv
import screed
import textwrap

from khmer.khmer_args import sanitize_help, KhmerArgumentParser


def get_parser():
    descr = "Display summary statistics for one or more FASTA/FASTQ files."
    epilog = """\
    Report number of bases, number of sequences, and average sequence length
    for one or more FASTA/FASTQ files; and report aggregate statistics at end.

    With :option:`-o`/:option:`--output`, the output will be saved to the
    specified file.

    Example::

        readstats.py tests/test-data/test-abund-read-2.fa
    """

    parser = KhmerArgumentParser(
        description=descr,
        epilog=textwrap.dedent(epilog))
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', dest='outfp', metavar="filename",
                        help="output file for statistics; defaults to stdout.",
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--csv', default=False, action='store_true',
                        help='Use the CSV format for the statistics, '
                        'including column headers.')
    return parser


class StatisticsOutput(object):  # pylint: disable=too-few-public-methods
    """
    Output statistics for several data files.

    The format of the output is determined by the formatter used.
    All statistics are aggregated and a summary is added to the data.
    """

    def __init__(self, formatter):
        self.formatter = formatter

    def __enter__(self):
        """Write header upon entry."""
        self.formatter.write_header()
        return self

    def append(self, basepairs, seqs, filename):
        """Append a new line for the given basepair num, sequences and file."""
        self.formatter.append(
            basepairs, seqs, basepairs / float(seqs), filename)

    def __exit__(self, exc_type, exc_value, traceback):
        """No exception? Finalize the formatter on exit."""
        if exc_type is None:
            self.formatter.finalize()


class CsvFormatter(object):
    """Format the statistics information as CSV."""

    headers = ['bp', 'seqs', 'avg_len', 'filename']

    def __init__(self, underlying_file):
        self.file = csv.writer(underlying_file)

    def write_header(self):
        """Add headers for the csv columns."""
        self.file.writerow(self.headers)

    def append(self, basepairs, seqs, avg_len, filename):
        """Append the data separated by comma."""
        self.file.writerow([basepairs, seqs, "%.1f" % avg_len, filename])

    def finalize(self):
        """No statistics since the CSV data is to be processed further."""
        pass


class StdFormatter(object):
    """Format the statistics in a human readable string."""

    def __init__(self, underlying_file):
        self.file = underlying_file
        self.bp_total = 0
        self.seqs_total = 0

    def write_header(self):
        """Write a header line."""
        self.file.write('---------------\n')

    def append(self, basepairs, seqs, avg_len, filename):
        """Append the data human readable."""
        self.bp_total += basepairs
        self.seqs_total += seqs
        self.file.write('%d bp / %d seqs; %.1f average length -- %s\n' %
                        (basepairs,
                         seqs,
                         avg_len,
                         filename))

    def finalize(self):
        """Add a summary with the accumulated data."""
        self.file.write('---------------\n')
        avg = self.bp_total / float(self.seqs_total)
        self.file.write('%d bp / %d seqs; %.1f average length -- total\n' %
                        (self.bp_total, self.seqs_total, avg))


def analyze_file(filename):
    """Run over the given file and count base pairs and sequences."""
    bps = 0
    seqs = 0
    input_iter = screed.open(filename)
    for record in input_iter:
        if seqs % 100000 == 0:
            print('...', filename, seqs, file=sys.stderr)
        bps += len(record.sequence)
        seqs += 1
    return bps, seqs


def main():
    """Main function - run when executed as a script."""
    args = sanitize_help(get_parser()).parse_args()

    statistics = []

    for filename in args.filenames:
        try:
            bps, seqs = analyze_file(filename)
        except (IOError, OSError, EOFError) as exc:
            print('ERROR in opening %s:' % filename, file=sys.stderr)
            print('     ', str(exc), file=sys.stderr)
            continue

        if seqs:
            statistics.append((bps, seqs, filename))
            avg = bps / float(seqs)
            msg = '%d bps / %d seqs; %.1f average length -- %s' % (bps,
                                                                   seqs,
                                                                   avg,
                                                                   filename)
            print('... found', msg, file=sys.stderr)
        else:
            print('No sequences found in %s' % filename, file=sys.stderr)

    if statistics:
        if args.csv:
            formatter = CsvFormatter(args.outfp)
        else:
            formatter = StdFormatter(args.outfp)
        with StatisticsOutput(formatter) as out:
            for stat in statistics:
                out.append(*stat)
    else:
        print('No sequences found in %d files' %
              len(args.filenames), file=args.outfp)


if __name__ == '__main__':
    main()
