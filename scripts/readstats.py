#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Display summary statistics for one or more FASTA/FASTQ files.

% scripts/readstats.py [ -o output.txt ] <file1> <file2>

Use '-h' for parameter help.
"""

import sys
import screed
import argparse
import textwrap


def get_parser():
    descr = "Display summary statistics for one or more FASTA/FASTQ files."
    epilog = ("""
    Report number of bases, number of sequences, and average sequence length
    for one or more FASTA/FASTQ files; and report aggregate statistics at end.

    With :option:`-o`/:options:`--output`, the output will be saved to the
    specified file.

    Example::

        readstats.py tests/test-data/test-abund-read-2.fa
    """)

    parser = argparse.ArgumentParser(description=descr,
                                     epilog=textwrap.dedent(epilog))
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', dest='outfp', metavar="filename",
                        help="output file for statistics; defaults to stdout.",
                        type=argparse.FileType('w'), default=None)

    return parser


def main():
    "Main function - run when executed as a script."

    parser = get_parser()
    args = parser.parse_args()

    total_bp = 0
    total_seqs = 0

    output = []
    for filename in args.filenames:
        bp = 0
        seqs = 0

        try:
            input_iter = screed.open(filename)
        except (IOError, OSError, EOFError) as exc:
            print >>sys.stderr, 'ERROR in opening %s:' % filename
            print >>sys.stderr, '     ', str(exc)
            continue

        for record in input_iter:
            if seqs % 100000 == 0:
                print >>sys.stderr, '...', filename, seqs
            bp += len(record.sequence)
            seqs += 1

        if seqs:
            avg_len = bp / float(seqs)
            s = '%d bp / %d seqs; %.1f average length -- %s' % (bp,
                                                                seqs,
                                                                avg_len,
                                                                filename)
            print >>sys.stderr, '... found', s
            output.append(s)

            total_bp += bp
            total_seqs += seqs
        else:
            print >>sys.stderr, 'No sequences found in %s' % filename

    if total_seqs:
        print >>args.outfp, '---------------'
        print >>args.outfp, "\n".join(output)
        print >>args.outfp, '---------------'
        print >>args.outfp, '%d bp / %d seqs; %.1f average length -- total' % \
            (total_bp, total_seqs, total_bp / float(total_seqs))
    else:
        print >>args.outfp, \
            'No sequences found in %d files' % len(args.filenames)


if __name__ == '__main__':
    main()
