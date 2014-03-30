#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Take an interleaved set of reads (/1 and /2), and extract them into separate
files (.1 and .2).

% scripts/split-paired-reads.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""
import screed
import sys
import os.path
import argparse
from khmer.file import check_file_status, check_space


def main():
    parser = argparse.ArgumentParser(
        description='Split interleaved reads into two files, left and right.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('infile')
    args = parser.parse_args()

    infile = args.infile

    check_file_status(infile)
    filenames = [infile]
    check_space(filenames)

    out1 = os.path.basename(infile) + '.1'
    out2 = os.path.basename(infile) + '.2'
    fp_out1 = open(out1, 'w')
    fp_out2 = open(out2, 'w')

    # is input file FASTQ or FASTA? Determine.
    is_fastq = False
    record = iter(screed.open(infile)).next()

    if hasattr(record, 'accuracy'):
        is_fastq = True

    counter1 = 0
    counter2 = 0
    for index, record in enumerate(screed.open(infile)):
        if index % 100000 == 0:
            print >> sys.stderr, '...', index

        name = record.name
        if name.endswith('/1'):
            if is_fastq:
                print >> fp_out1, '@%s\n%s\n+\n%s' % (record.name,
                                                      record.sequence,
                                                      record.accuracy)
            else:
                print >> fp_out1, '>%s\n%s' % (record.name, record.sequence,)
            counter1 += 1
        elif name.endswith('/2'):
            if is_fastq:
                print >> fp_out2, '@%s\n%s\n+\n%s' % (record.name,
                                                      record.sequence,
                                                      record.accuracy)
            else:
                print >> fp_out2, '>%s\n%s' % (record.name, record.sequence,)
            counter2 += 1

    print >> sys.stderr, "DONE; split %d sequences (%d left, %d right)" % \
        (index + 1, counter1, counter2)
    print >> sys.stderr, "/1 reads in %s" % out1
    print >> sys.stderr, "/2 reads in %s" % out2

if __name__ == '__main__':
    main()
