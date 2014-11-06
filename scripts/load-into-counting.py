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

import os
import sys
import threading
import textwrap
import khmer
from khmer.khmer_args import build_counting_args, report_on_config, info,\
    add_threading_args
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

        load-into-counting.py -k 20 -x 5e7 out.kh data/100k-filtered.fa

    Multiple threads can be used to accelerate the process, if you have extra
    cores to spare.

    Example::

        load-into-counting.py -k 20 -x 5e7 -T 4 out.kh data/100k-filtered.fa
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
    parser.add_argument('--machine-readable-info', '-m', action='store_true',
                        help="Should we also create a .info.tsv file "
                        "containing a machine readable run summary?")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
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

    print >>sys.stderr, 'Saving k-mer counting table',  base
    print >>sys.stderr, 'Loading kmers from sequences in', repr(filenames)

    # clobber the '.info' file now, as we always open in append mode below
    if os.path.exists(base + '.info'):
        os.remove(base + '.info')

    print >>sys.stderr, 'making k-mer counting table'
    htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                     args.n_tables, args.threads)
    htable.set_use_bigcount(args.bigcount)

    config = khmer.get_config()
    config.set_reads_input_buffer_size(args.threads * 64 * 1024)

    filename = None

    for index, filename in enumerate(filenames):

        rparser = khmer.ReadParser(filename, args.threads)
        threads = []
        print >>sys.stderr, 'consuming input', filename
        for _ in xrange(args.threads):
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
            print >>sys.stderr, 'mid-save', base
            htable.save(base)
        with open(base + '.info', 'a') as info_fh:
            print >> info_fh, 'through', filename

    n_kmers = htable.n_occupied()
    if args.report_total_kmers:
        print >> sys.stderr, 'Total number of k-mers:', n_kmers
        with open(base + '.info', 'a') as info_fp:
            print >>info_fp, 'total distinct k-mers:', n_kmers

    print >>sys.stderr, 'saving', base
    htable.save(base)

    fp_rate = khmer.calc_expected_collisions(htable)

    with open(base + '.info', 'a') as info_fp:
        print >> sys.stderr, "Writing run information to", base + '.info'
        print >> info_fp, 'fp rate estimated to be %1.3f\n' % fp_rate


    if args.machine_readable_info:
        tsv_file = base + '.info.tsv'
        print >> sys.stderr, "Writing machine-readable stats to", tsv_file
        if os.path.exists(tsv_file):
            os.remove(tsv_file)
        with open(tsv_file, 'a') as tsv_fh:
            tsv_fh.write("name\tfpr\tkmers\tfiles\n")
            tsv_fh.write("{base:s}\t{fpr:1.3f}\t{kmers:d}\t{files:s}\n".format(
                base=os.path.basename(base), fpr=fp_rate, kmers=n_kmers,
                files=";".join(filenames)))

    print >> sys.stderr, 'fp rate estimated to be %1.3f' % fp_rate

    # Change 0.2 only if you really grok it.  HINT: You don't.
    if fp_rate > 0.20:
        print >> sys.stderr, "**"
        print >> sys.stderr, "** ERROR: the k-mer counting table is too small",
        print >> sys.stderr, "for this data set. Increase tablesize/# tables."
        print >> sys.stderr, "**"
        sys.exit(1)

    print >>sys.stderr, 'DONE.'
    print >>sys.stderr, 'wrote to:', base + '.info'

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
