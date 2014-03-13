#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Count the median/avg k-mer abundance for each sequence in the input file,
based on the k-mer counts in the given counting hash.  Can be used to
estimate expression levels (mRNAseq) or coverage (genomic/metagenomic).

% scripts/count-median.py <htname> <input seqs> <output counts>

Use '-h' for parameter help.

The output file contains sequence id, median, average, stddev, and seq length.

NOTE: All 'N's in the input sequences are converted to 'G's.
"""
import screed
import khmer
import argparse
from khmer.file_api import check_file_status, check_space
#


def main():
    parser = argparse.ArgumentParser(
        description='Count k-mers summary stats for sequences')

    parser.add_argument('htfile')
    parser.add_argument('input')
    parser.add_argument('output')

    args = parser.parse_args()

    htfile = args.htfile
    input_filename = args.input
    output_filename = args.output

    infiles = [htfile, input_filename]
    for infile in infiles:
        check_file_status(infile)

    check_space(infiles)

    print 'loading counting hash from', htfile
    htable = khmer.load_counting_hash(htfile)
    ksize = htable.ksize()

    print 'writing to', output_filename
    output = open(output_filename, 'w')

    for record in screed.open(input_filename):
        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'G')

        if ksize <= len(seq):
            medn, ave, stdev = htable.get_median_count(seq)
            print >> output, record.name, medn, ave, stdev, len(seq)

if __name__ == '__main__':
    main()
