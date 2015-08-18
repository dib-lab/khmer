#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,invalid-name
"""
Sequence trimming by abundance w/o countgraph.

Trim sequences at k-mers of the given abundance for the given file,
without loading a prebuilt countgraph.  Output sequences will be
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
from khmer import khmer_args
from khmer.khmer_args import (build_counting_args, report_on_config,
                              add_threading_args, info, calculate_graphsize,
                              sanitize_help)
from khmer.kfile import (check_input_files, check_space,
                         check_space_for_graph,
                         add_output_compression_type,
                         get_file_writer)

DEFAULT_CUTOFF = 2


def get_parser():
    epilog = """\
    Trimmed sequences will be placed in
    ``${input_sequence_filename}.abundfilt``.

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
    parser.add_argument('--savegraph', metavar="filename", default='',
                        help="If present, the name of the file to save the "
                        "k-mer countgraph to")
    parser.add_argument('datafile', metavar='input_sequence_filename',
                        help="FAST[AQ] sequence file to trim")
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    add_output_compression_type(parser)
    return parser


def main():
    info('filter-abund-single.py', ['counting', 'SeqAn'])
    args = sanitize_help(get_parser()).parse_args()

    check_input_files(args.datafile, args.force)
    check_space([args.datafile], args.force)

    if args.savegraph:
        tablesize = calculate_graphsize(args, 'countgraph')
        check_space_for_graph(args.savegraph, tablesize, args.force)

    report_on_config(args)

    print('making countgraph', file=sys.stderr)
    graph = khmer_args.create_countgraph(args)

    # first, load reads into graph
    rparser = khmer.ReadParser(args.datafile)
    threads = []
    print('consuming input, round 1 --', args.datafile, file=sys.stderr)
    for _ in range(args.threads):
        cur_thread = \
            threading.Thread(
                target=graph.consume_fasta_with_reads_parser,
                args=(rparser, )
            )
        threads.append(cur_thread)
        cur_thread.start()

    for _ in threads:
        _.join()

    print('Total number of unique k-mers: {0}'.format(
        graph.n_unique_kmers()), file=sys.stderr)

    fp_rate = khmer.calc_expected_collisions(graph, args.force)
    print('fp rate estimated to be %1.3f' % fp_rate, file=sys.stderr)

    # now, trim.

    # the filtering function.
    def process_fn(record):
        name = record.name
        seq = record.sequence
        seqN = seq.replace('N', 'A')

        _, trim_at = graph.trim_on_abundance(seqN, args.cutoff)

        if trim_at >= args.ksize:
            # be sure to not to change the 'N's in the trimmed sequence -
            # so, return 'seq' and not 'seqN'.
            return name, seq[:trim_at]

        return None, None

    # the filtering loop
    print('filtering', args.datafile, file=sys.stderr)
    outfile = os.path.basename(args.datafile) + '.abundfilt'
    outfile = open(outfile, 'wb')
    outfp = get_file_writer(outfile, args.gzip, args.bzip)

    tsp = ThreadedSequenceProcessor(process_fn)
    tsp.start(verbose_loader(args.datafile), outfp)

    print('output in', outfile, file=sys.stderr)

    if args.savegraph:
        print('Saving k-mer countgraph filename',
              args.savegraph, file=sys.stderr)
        print('...saving to', args.savegraph, file=sys.stderr)
        graph.save(args.savegraph)
    print('wrote to: ', outfile, file=sys.stderr)

if __name__ == '__main__':
    main()
