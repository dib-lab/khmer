#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Produce the k-mer abundance distribution for the given file.

% python scripts/abundance-dist.py [ -z -s ] <htname> <data> <histout>

Use '-h' for parameter help.
"""
import sys
import khmer
import argparse
import os

#  Import fileapi from sandbox - temporary arrangement
current_file_path = os.path.realpath(__file__)
current_folder = os.path.dirname(current_file_path)
parent_folder = os.path.dirname(current_folder)
sandbox_folder = os.path.join(parent_folder, 'sandbox')
sys.path.append(sandbox_folder)

import fileApi

def main():
    parser = argparse.ArgumentParser(
        description="Output k-mer abundance distribution.")

    parser.add_argument('hashname')
    parser.add_argument('datafile')
    parser.add_argument('histout')

    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')

    args = parser.parse_args()
    hashfile = args.hashname
    datafile = args.datafile
    histout = args.histout
    
    # Check if input files exist
    infiles = [hashfile, datafile]
    for infile in infiles:
        fileApi.check_file_status(infile)
    
    # Check free space
    fileApi.check_space(infiles)

    print 'hashtable from', hashfile
    ht = khmer.load_counting_hash(hashfile)

    K = ht.ksize()
    sizes = ht.hashsizes()
    tracking = khmer._new_hashbits(K, sizes)

    print 'K:', K
    print 'HT sizes:', sizes
    print 'outputting to', histout

    if os.path.exists(histout):
        if not args.squash_output:
            print >>sys.stderr, 'ERROR: %s exists; not squashing.' % histout
            sys.exit(-1)

        print '** squashing existing file %s' % histout

    print 'preparing hist...'
    z = ht.abundance_distribution(datafile, tracking)
    total = sum(z)

    if 0 == total:
        print >>sys.stderr, \
            "ERROR: abundance distribution is uniformly zero; " \
            "nothing to report."
        print >>sys.stderr, "\tPlease verify that the input files are valid."
        sys.exit(-1)

    fp = open(histout, 'w')

    sofar = 0
    for n, i in enumerate(z):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print >>fp, n, i, sofar, round(frac, 3)

        if sofar == total:
            break

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
