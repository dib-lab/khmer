#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
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
from khmer.khmer_logger import (configure_logging, log_info, log_error,
                                log_warn)


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
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
    add_output_compression_type(parser)
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()
    if not args.quiet:
        info('filter-abund-single.py', ['counting', 'SeqAn'])

    configure_logging(args.quiet)
    check_input_files(args.datafile, args.force)
    check_space([args.datafile], args.force)

    if args.savegraph:
        tablesize = calculate_graphsize(args, 'countgraph')
        check_space_for_graph(args.savegraph, tablesize, args.force)

    report_on_config(args)

    log_info('making countgraph')
    graph = khmer_args.create_countgraph(args)

    # first, load reads into graph
    rparser = khmer.ReadParser(args.datafile)
    threads = []
    log_info('consuming input, round 1 -- {datafile}', datafile=args.datafile)
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

    log_info('Total number of unique k-mers: {nk}', nk=graph.n_unique_kmers())

    fp_rate = khmer.calc_expected_collisions(graph, args.force)
    log_info('fp rate estimated to be {fpr:1.3f}', fpr=fp_rate)

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
    log_info('filtering {datafile}', datafile=args.datafile)
    outfile = os.path.basename(args.datafile) + '.abundfilt'
    outfile = open(outfile, 'wb')
    outfp = get_file_writer(outfile, args.gzip, args.bzip)

    tsp = ThreadedSequenceProcessor(process_fn, verbose=not args.quiet)
    tsp.start(verbose_loader(args.datafile), outfp)

    log_info('output in {outfile}', outfile=outfile.name)

    if args.savegraph:
        log_info('Saving k-mer countgraph filename {graph}',
                 graph=args.savegraph)
        graph.save(args.savegraph)

if __name__ == '__main__':
    main()
