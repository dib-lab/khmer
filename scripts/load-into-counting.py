#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name
"""
Build a counting Bloom filter from the given sequences, save in <htname>.

% load-into-counting.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import threading
import textwrap
import khmer
from khmer.khmer_args import build_counting_args, report_on_config, info
from khmer.threading_args import add_threading_args
from khmer.file import check_file_status, check_space
from khmer.file import check_space_for_hashtable


def get_parser():
    epilog = """
    Note: with :option:`-b` the output will be the exact size of the
    k-mer counting table and this script will use a constant amount of memory.
    In exchange k-mer counts will stop at 255. The memory usage of this script
    with :option:`-b` will be about 1.15x the product of the :option:`-x` and
    :option:`-N` numbers.

    Example::

        load_into_counting.py -k 20 -x 5e7 out.kh data/100k-filtered.fa

    Multiple threads can be used to accelerate the process, if you have extra
    cores to spare.

    Example::

        load_into_counting.py -k 20 -x 5e7 -T 4 out.kh data/100k-filtered.fa
    """

    parser = build_counting_args("Build a k-mer counting table from the given"
                                 " sequences.", epilog=textwrap.dedent(epilog))
    add_threading_args(parser)
    parser.add_argument('output_countingtable_filename', help="The name of the"
                        " file to write the k-mer counting table to.")
    parser.add_argument('input_sequence_filename', nargs='+',
                        help="The names of one or more FAST[AQ] input "
                        "sequence files.")
    parser.add_argument('-b', '--no-bigcount', dest='bigcount', default=True,
                        action='store_false',
                        help='Do not count k-mers past 255')
    return parser


def main():

    info('load-into-counting.py', ['counting'])
    args = get_parser().parse_args()
    report_on_config(args)

    base = args.output_countingtable_filename
    filenames = args.input_sequence_filename

    for name in args.input_sequence_filename:
        check_file_status(name)

    check_space(args.input_sequence_filename)
    check_space_for_hashtable(args.n_tables * args.min_tablesize)

    print 'Saving k-mer counting table to %s' % base
    print 'Loading kmers from sequences in %s' % repr(filenames)

    print 'making k-mer counting table'
    htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                     args.n_tables, args.n_threads)
    htable.set_use_bigcount(args.bigcount)

    config = khmer.get_config()
    config.set_reads_input_buffer_size(args.n_threads * 64 * 1024)

    for index, filename in enumerate(filenames):

        rparser = khmer.ReadParser(filename, args.n_threads)
        threads = []
        print 'consuming input', filename
        for _ in xrange(args.n_threads):
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
            check_space_for_hashtable(args.n_tables * args.min_tablesize)
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
        print >> sys.stderr, ("** ERROR: the k-mer counting table is too small"
                              " this data set.  Increase tablesize/# tables.")
        print >> sys.stderr, "**"
        sys.exit(1)

    print 'DONE.'

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
