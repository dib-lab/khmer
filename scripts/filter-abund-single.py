#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Trim sequences at k-mers of the given abundance for the given file,
without loading a prebuilt counting hash.  Output sequences will be
placed in 'infile.abundfilt'.

% python scripts/filter-abund-single.py <data>

Use '-h' for parameter help.
"""
import os
import khmer
import threading
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer.counting_args import build_construct_args, report_on_config
from khmer.threading_args import add_threading_args

from khmer.file import check_file_status, check_space
from khmer.file import check_space_for_hashtable
#

DEFAULT_CUTOFF = 2


def main():
    parser = build_construct_args(
        "Filter k-mers at the given abundance (inmem version).")
    add_threading_args(parser)

    parser.add_argument('--cutoff', '-C', dest='cutoff',
                        default=DEFAULT_CUTOFF, type=int,
                        help="Trim at k-mers below this abundance.")
    parser.add_argument('--savehash', dest='savehash', default='')
    parser.add_argument('datafile')

    args = parser.parse_args()
    report_on_config(args)

    ksize = args.ksize
    min_hashsize = args.min_hashsize
    n_hashes = args.n_hashes
    n_threads = int(args.n_threads)
    filename = args.datafile

    check_file_status(filename)
    infiles = [filename]
    check_space(infiles)
    if args.savehash:
        check_space_for_hashtable(ksize * min_hashsize)

    config = khmer.get_config()
    bufsz = config.get_reads_input_buffer_size()
    config.set_reads_input_buffer_size(n_threads * 64 * 1024)

    print 'making hashtable'
    htable = khmer.new_counting_hash(ksize, min_hashsize, n_hashes, n_threads)

    # first, load reads into hash table
    rparser = khmer.ReadParser(filename, n_threads)
    threads = []
    print 'consuming input, round 1 --', filename
    for tnum in xrange(n_threads):
        cur_thread = \
            threading.Thread(
                target=htable.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(cur_thread)
        cur_thread.start()

    for _ in threads:
        _.join()

    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be %1.3f' % fp_rate

    # now, trim.

    # the filtering function.
    def process_fn(record):
        name = record['name']
        seq = record['sequence']
        if 'N' in seq:
            return None, None

        trim_seq, trim_at = htable.trim_on_abundance(seq, args.cutoff)

        if trim_at >= ksize:
            return name, trim_seq

        return None, None

    # the filtering loop
    print 'filtering', filename
    outfile = os.path.basename(filename) + '.abundfilt'
    outfp = open(outfile, 'w')

    tsp = ThreadedSequenceProcessor(process_fn)
    tsp.start(verbose_loader(filename), outfp)

    print 'output in', outfile

    if args.savehash:
        print 'Saving hashfile', args.savehash
        print '...saving to', args.savehash
        htable.save(args.savehash)

if __name__ == '__main__':
    main()
