#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Sequence trimming by abundance w/o counting table.

Trim sequences at k-mers of the given abundance for the given file,
without loading a prebuilt counting table.  Output sequences will be
placed in 'infile.abundfilt'.

% python scripts/filter-abund-single.py <data>

Use '-h' for parameter help.
"""
from __future__ import print_function
import os
import sys
import khmer
import threading
import textwrap
from khmer.thread_utils import ThreadedSequenceProcessor, verbose_loader
from khmer.khmer_args import (build_counting_args, report_on_config,
                              add_threading_args, info)
from khmer.kfile import (check_input_files, check_space,
                         check_space_for_hashtable)
#
DEFAULT_CUTOFF = 2


def get_parser():
    epilog = """
    Trimmed sequences will be placed in ${input_sequence_filename}.abundfilt.

    This script is constant memory.

    To trim reads based on k-mer abundance across multiple files, use
    :program:`load-into-counting.py` and :program:`filter-abund.py`.

    Example::

        filter-abund-single.py -k 20 -x 5e7 -C 2 data/100k-filtered.fa
    """
    parser = build_counting_args(
        descr="Trims sequences at a minimum k-mer abundance "
        "(in memory version).", epilog=textwrap.dedent(epilog))
    add_threading_args(parser)

    parser.add_argument('--cutoff', '-C', default=DEFAULT_CUTOFF, type=int,
                        help="Trim at k-mers below this abundance.")
    parser.add_argument('--savetable', metavar="filename", default='',
                        help="If present, the name of the file to save the "
                        "k-mer counting table to")
    parser.add_argument('datafile', metavar='input_sequence_filename',
                        help="FAST[AQ] sequence file to trim")
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers to stderr")
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    return parser


def main():
    info('filter-abund-single.py', ['counting', 'SeqAn'])
    args = get_parser().parse_args()
    check_input_files(args.datafile, args.force)
    check_space([args.datafile], args.force)
    if args.savetable:
        check_space_for_hashtable(
            args.n_tables * args.min_tablesize, args.force)
    report_on_config(args)

    print('making k-mer counting table', file=sys.stderr)
    htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                     args.n_tables)

    # first, load reads into hash table
    rparser = khmer.ReadParser(args.datafile)
    threads = []
    print('consuming input, round 1 --', args.datafile, file=sys.stderr)
    for _ in range(args.threads):
        cur_thread = \
            threading.Thread(
                target=htable.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(cur_thread)
        cur_thread.start()

    for _ in threads:
        _.join()

    if args.report_total_kmers:
        print('Total number of unique k-mers: {0}'.format(
            htable.n_unique_kmers()), file=sys.stderr)

    fp_rate = khmer.calc_expected_collisions(htable, args.force)
    print('fp rate estimated to be %1.3f' % fp_rate, file=sys.stderr)

    # now, trim.

    # the filtering function.
    def process_fn(record):
        name = record.name
        seq = record.sequence
        seqN = seq.replace('N', 'A')

        _, trim_at = htable.trim_on_abundance(seqN, args.cutoff)

        if trim_at >= args.ksize:
            # be sure to not to change the 'N's in the trimmed sequence -
            # so, return 'seq' and not 'seqN'.
            return name, seq[:trim_at]

        return None, None

    # the filtering loop
    print('filtering', args.datafile, file=sys.stderr)
    outfile = os.path.basename(args.datafile) + '.abundfilt'
    outfp = open(outfile, 'w')

    tsp = ThreadedSequenceProcessor(process_fn)
    tsp.start(verbose_loader(args.datafile), outfp)

    print('output in', outfile, file=sys.stderr)

    if args.savetable:
        print('Saving k-mer counting table filename',
              args.savetable, file=sys.stderr)
        print('...saving to', args.savetable, file=sys.stderr)
        htable.save(args.savetable)
    print('wrote to: ', outfile, file=sys.stderr)

if __name__ == '__main__':
    main()
