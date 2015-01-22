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


def main():
    descr = "Display summary statistics for one or more FASTA/FASTQ files."
    epilog = """ """
    
    parser = argparse.ArgumentParser(description=descr,
                                     epilog=textwrap.dedent(epilog))
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', dest='outfp',
                        help="output file for statistics; defaults to stdout.",
                        type=argparse.FileType('w'), default=None)
    args = parser.parse_args()

    total_bp = 0
    total_seqs = 0

    output = []
    for filename in args.filenames:
        bp = 0
        seqs = 0
        for record in screed.open(filename):
            if seqs % 100000 == 0:
                print >>sys.stderr, '...', filename, seqs
            bp += len(record.sequence)
            seqs += 1

        if seqs == 0:
            print >>sys.stderr, 'No sequences found in %s' % filename
        else:
            avg_len = bp / float(seqs)
            s = '%d bp / %d seqs; %.1f average length -- %s' % (bp,
                                                                seqs,
                                                                avg_len,
                                                                filename)
            print >>sys.stderr, '... found', s
            output.append(s)

            total_bp += bp
            total_seqs += seqs

    if total_seqs == 0:
        print >>args.outfp, \
            'No sequences found in %d files' % len(args.filenames)
    else:
        print >>args.outfp, '---------------'
        print >>args.outfp, "\n".join(output)
        print >>args.outfp, '---------------'
        print >>args.outfp, '%d bp / %d seqs; %.1f average length -- total' % \
            (total_bp, total_seqs, total_bp / float(total_seqs))


if __name__ == '__main__':
    main()
