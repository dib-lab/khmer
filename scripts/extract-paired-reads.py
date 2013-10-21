#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Take a file containing a mixture of interleaved and orphaned reads, and
extract them into separate files (.pe and .se).

% scripts/extract-paired-reads.py <infile>

Reads FASTQ and FASTA input, retains format for output.
"""
import screed
import sys
import os.path
import argparse


def is_pair(name1, name2):
    if name1.endswith('/1') and name2.endswith('/2'):
        s1 = name1.split('/')[0]
        s2 = name2.split('/')[0]
        if s1 == s2:
            assert(s1)
            return True

    return False


def output_pair(r1, r2):
    if hasattr(r1, 'accuracy'):
        return "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n" % \
            (r1.name, r1.sequence, r1.accuracy,
             r2.name, r2.sequence, r2.accuracy)
    else:
        return ">%s\n%s\n>%s\n%s\n" % (r1.name, r1.sequence, r2.name,
                                       r2.sequence)


def output_single(r):
    if hasattr(r, 'accuracy'):
        return "@%s\n%s\n+\n%s\n" % (r.name, r.sequence, r.accuracy)
    else:
        return ">%s\n%s\n" % (r.name, r.sequence)


def main():
    parser = argparse.ArgumentParser(
        description=
        'Take a mixture of reads and split into pairs and orphans.')

    parser.add_argument('infile')
    args = parser.parse_args()

    infile = args.infile
    outfile = os.path.basename(infile)
    if len(sys.argv) > 2:
        outfile = sys.argv[2]

    single_fp = open(outfile + '.se', 'w')
    paired_fp = open(outfile + '.pe', 'w')

    print 'reading file "%s"' % infile
    print 'outputting interleaved pairs to "%s.pe"' % outfile
    print 'outputting orphans to "%s.se"' % outfile

    last_record = None
    last_name = None

    n_pe = 0
    n_se = 0

    record = None
    for n, record in enumerate(screed.open(sys.argv[1])):
        if n % 100000 == 0 and n > 0:
            print '...', n
        name = record['name'].split()[0]

        if last_record:
            if is_pair(last_name, name):
                paired_fp.write(output_pair(last_record, record))
                name, record = None, None
                n_pe += 1
            else:
                single_fp.write(output_single(last_record))
                n_se += 1

        last_name = name
        last_record = record

    if last_record:
        if is_pair(last_name, name):
            paired_fp.write(output_pair(last_record, record))
            name, record = None, None
            n_pe += 1
        else:
            single_fp.write(output_single(last_record))
            name, record = None, None
            n_se += 1

    if record:
        single_fp.write(output_single(record))
        n_se += 1

    single_fp.close()
    paired_fp.close()

    if n_pe == 0:
        raise Exception("no paired reads!? check file formats...")

    print 'DONE; read %d sequences, %d pairs and %d singletons' % \
          (n + 1, n_pe, n_se)

if __name__ == '__main__':
    main()
