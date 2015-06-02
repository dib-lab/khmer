#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name
"""
Build a counting Bloom filter from the given sequences, save in <htname>.

% load-into-counting.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
from __future__ import print_function, unicode_literals

import json
import os
import sys
import threading
import textwrap
import khmer
from khmer.khmer_args import build_counting_args, report_on_config, info,\
    add_threading_args
from khmer.kfile import check_file_writable
from khmer.kfile import check_input_files, check_space
from khmer.kfile import check_space_for_hashtable


def get_parser():
    epilog = """
    Note: with :option:`-b` the output will be the exact size of the
    k-mer counting table and this script will use a constant amount of memory.
    In exchange k-mer counts will stop at 255. The memory usage of this script
    with :option:`-b` will be about 1.15x the product of the :option:`-x` and
    :option:`-N` numbers.

    Example::

        load-into-counting.py -k 20 -x 5e7 out.ct data/100k-filtered.fa

    Multiple threads can be used to accelerate the process, if you have extra
    cores to spare.

    Example::

        load-into-counting.py -k 20 -x 5e7 -T 4 out.ct data/100k-filtered.fa
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
                        action='store_false', help="The default behaviour is "
                        "to count past 255 using bigcount. This flag turns "
                        "bigcount off, limiting counts to 255.")
    parser.add_argument('--summary-info', '-s', type=str, default=None,
                        metavar="FORMAT", choices=[str('json'), str('tsv')],
                        help="What format should the machine readable run "
                        "summary be in? (json or tsv, disabled by default)")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():

    info('load-into-counting.py', ['counting', 'SeqAn'])

    args = get_parser().parse_args()
    report_on_config(args)

    base = args.output_countingtable_filename
    filenames = args.input_sequence_filename

    for name in args.input_sequence_filename:
        check_input_files(name, args.force)

    check_space(args.input_sequence_filename, args.force)
    check_space_for_hashtable(args.n_tables * args.min_tablesize, args.force)

    check_file_writable(base)
    check_file_writable(base + ".info")

    print('Saving k-mer counting table to %s' % base, file=sys.stderr)
    print('Loading kmers from sequences in %s' %
          repr(filenames), file=sys.stderr)

    # clobber the '.info' file now, as we always open in append mode below
    if os.path.exists(base + '.info'):
        os.remove(base + '.info')

    print('making k-mer counting table', file=sys.stderr)
    htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                     args.n_tables)
    htable.set_use_bigcount(args.bigcount)

    filename = None

    total_num_reads = 0

    for index, filename in enumerate(filenames):

        rparser = khmer.ReadParser(filename)
        threads = []
        print('consuming input', filename, file=sys.stderr)
        for _ in range(args.threads):
            cur_thrd = \
                threading.Thread(
                    target=htable.consume_fasta_with_reads_parser,
                    args=(rparser, )
                )
            threads.append(cur_thrd)
            cur_thrd.start()

        for thread in threads:
            thread.join()

        if index > 0 and index % 10 == 0:
            check_space_for_hashtable(args.n_tables * args.min_tablesize,
                                      args.force)
            print('mid-save', base, file=sys.stderr)
            htable.save(base)
        with open(base + '.info', 'a') as info_fh:
            print('through', filename, file=info_fh)
        total_num_reads += rparser.num_reads

    n_kmers = htable.n_unique_kmers()
    if args.report_total_kmers:
        print('Total number of unique k-mers:', n_kmers, file=sys.stderr)
        with open(base + '.info', 'a') as info_fp:
            print('Total number of unique k-mers:', n_kmers, file=info_fp)

    print('saving', base, file=sys.stderr)
    htable.save(base)

    # Change max_false_pos=0.2 only if you really grok it. HINT: You don't
    fp_rate = \
        khmer.calc_expected_collisions(htable, args.force, max_false_pos=.2)

    with open(base + '.info', 'a') as info_fp:
        print('fp rate estimated to be %1.3f\n' % fp_rate, file=info_fp)

    if args.summary_info:
        mr_fmt = args.summary_info.lower()
        mr_file = base + '.info.' + mr_fmt
        print("Writing summmary info to", mr_file, file=sys.stderr)
        with open(mr_file, 'w') as mr_fh:
            if mr_fmt == 'json':
                mr_data = {
                    "ht_name": os.path.basename(base),
                    "fpr": fp_rate,
                    "num_kmers": n_kmers,
                    "files": filenames,
                    "mrinfo_version": "0.2.0",
                    "num_reads": total_num_reads,
                }
                json.dump(mr_data, mr_fh)
                mr_fh.write('\n')
            elif mr_fmt == 'tsv':
                mr_fh.write("ht_name\tfpr\tnum_kmers\tnum_reads\tfiles\n")
                vals = [
                    os.path.basename(base),
                    "{:1.3f}".format(fp_rate),
                    str(n_kmers),
                    str(total_num_reads),
                    ";".join(filenames),
                ]
                mr_fh.write("\t".join(vals) + "\n")

    print('fp rate estimated to be %1.3f' % fp_rate, file=sys.stderr)

    print('DONE.', file=sys.stderr)
    print('wrote to:', base + '.info', file=sys.stderr)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
