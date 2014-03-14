#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Build a counting Bloom filter from the given sequences, save in <htname>.

% python scripts/load-into-counting.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import threading
import khmer
from khmer.khmer_args import build_counting_args, report_on_config
from khmer.threading_args import add_threading_args
from khmer.file_api import check_file_status, check_space
from khmer.file_api import check_space_for_hashtable
#


def main():
    parser = build_counting_args()
    add_threading_args(parser)
    parser.add_argument('output_filename')
    parser.add_argument('input_filenames', nargs='+')
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')

    args = parser.parse_args()
    report_on_config(args)

    ksize = args.ksize
    min_hashsize = args.min_hashsize
    n_hashes = args.n_hashes

    base = args.output_filename
    filenames = args.input_filenames
    n_threads = int(args.n_threads)

    for _ in args.input_filenames:
        check_file_status(_)

    check_space(args.input_filenames)
    check_space_for_hashtable(args.ksize * args.min_hashsize)

    print 'Saving hashtable to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)

    #

    print 'making hashtable'
    htable = khmer.new_counting_hash(ksize, min_hashsize, n_hashes, n_threads)
    htable.set_use_bigcount(args.bigcount)

    config = khmer.get_config()
    bufsz = config.get_reads_input_buffer_size()
    config.set_reads_input_buffer_size(n_threads * 64 * 1024)

    for index, filename in enumerate(filenames):

        rparser = khmer.ReadParser(filename, n_threads)
        threads = []
        print 'consuming input', filename
        for tnum in xrange(n_threads):
            cur_thrd = \
                threading.Thread(
                    target=htable.consume_fasta_with_reads_parser,
                    args=(rparser, )
                )
            threads.append(cur_thrd)
            cur_thrd.start()

        for _ in threads:
            _.join()

        if index > 0 and index % 10 == 0:
            check_space_for_hashtable(ksize * min_hashsize)
            print 'mid-save', base
            htable.save(base)
            open(base + '.info', 'w').write('through %s' % filename)

    print 'saving', base
    htable.save(base)

    info_fp = open(base + '.info', 'w')
    info_fp.write('through end: %s\n' % filename)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be %1.3f' % fp_rate
    print >> info_fp, 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >> sys.stderr, "**"
        print >> sys.stderr, "** ERROR: the counting hash is too small for"
        print >> sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >> sys.stderr, "**"
        sys.exit(1)

    print 'DONE.'

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
