#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Produce the k-mer abundance distribution for the given file, without
loading a prebuilt counting hash.

% python scripts/abundance-dist-single.py <data> <histout>

Use '-h' for parameter help.
"""
import sys
import khmer
import threading
from khmer.counting_args import build_construct_args, report_on_config
from khmer.threading_args import add_threading_args

# Add sandbox to path - when fileApi is moved to 
# scripts/, this can be removed
current_file_path = os.path.realpath(__file__)
current_folder = os.path.dirname(current_file_path)
parent_folder = os.path.dirname(current_folder)
sandbox_folder = os.path.join(parent_folder, 'sandbox')
sys.path.append(sandbox_folder)

import fileApi
import datetime
import time

def main():
    parser = build_construct_args(
        "Output k-mer abundance distribution (single file version).")
    add_threading_args(parser)

    parser.add_argument('datafile')
    parser.add_argument('histout')

    parser.add_argument('-z', '--no-zero', dest='output_zero', default=True,
                        action='store_false',
                        help='Do not output 0-count bins')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    parser.add_argument('-s', '--squash', dest='squash_output', default=False,
                        action='store_true',
                        help='Overwrite output file if it exists')
    parser.add_argument('--savehash', dest='savehash', default='')

    args = parser.parse_args()
    report_on_config(args)

    print 'making hashtable'
    counting_hash = khmer.new_counting_hash(args.ksize, args.min_hashsize,
                                            args.n_hashes,
                                            args.n_threads)
    counting_hash.set_use_bigcount(args.bigcount)

    print 'building tracking ht'
    tracking = khmer.new_hashbits(counting_hash.ksize(), args.min_hashsize,
                                  args.n_hashes)

    print 'kmer_size:', counting_hash.ksize()
    print 'counting hash sizes:', counting_hash.hashsizes()
    print 'outputting to', args.histout

    khmer.get_config().set_reads_input_buffer_size(args.n_threads * 64 * 1024)

    # start loading
    rparser = khmer.ReadParser(args.datafile, args.n_threads)
    threads = []
    print 'consuming input, round 1 --', args.datafile
    for _ in xrange(args.n_threads):
        thread = \
            threading.Thread(
                target=counting_hash.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    abundance_lists = []

    def __do_abundance_dist__(read_parser):
        abundances = counting_hash.abundance_distribution_with_reads_parser(
            read_parser, tracking)
        abundance_lists.append(abundances)

    print 'preparing hist from %s...' % args.datafile
    rparser = khmer.ReadParser(args.datafile, args.n_threads)
    threads = []
    print 'consuming input, round 2 --', args.datafile
    for _ in xrange(args.n_threads):
        thread = \
            threading.Thread(
                target=__do_abundance_dist__,
                args=(rparser,)
            )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    assert len(abundance_lists) == args.n_threads, len(abundance_lists)
    abundance = {}
    for abundance_list in abundance_lists:
        for i, count in enumerate(abundance_list):
            abundance[i] = abundance.get(i, 0) + count

    total = sum(abundance.values())

    if 0 == total:
        print >>sys.stderr, \
            "ERROR: abundance distribution is uniformly zero; " \
            "nothing to report."
        print >>sys.stderr, "\tPlease verify that the input files are valid."
        sys.exit(1)

    hist_fp = open(args.histout, 'w')

    sofar = 0
    for _, i in sorted(abundance.items()):
        if i == 0 and not args.output_zero:
            continue

        sofar += i
        frac = sofar / float(total)

        print >>hist_fp, _, i, sofar, round(frac, 3)

        if sofar == total:
            break

    if args.savehash:
        print 'Saving hashfile', args.savehash
        print '...saving to', args.savehash
        counting_hash.save(args.savehash)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
