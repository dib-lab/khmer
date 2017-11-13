#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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

# pylint: disable=invalid-name,missing-docstring
"""
Count saturation curve for reads with a coverage of 1, but collect
reads whether or not they have high coverage.  This is better for
assessing saturation of (esp) low-coverage data sets.
"""

import sys
import screed
import os
import khmer
import textwrap

from khmer.khmer_args import (build_counting_args, add_loadgraph_args,
                              report_on_config, info, create_countgraph,
                              sanitize_help)
import argparse
from khmer.kfile import (check_space, check_space_for_graph,
                         check_valid_file_exists)
DEFAULT_DESIRED_COVERAGE = 1


# Iterate a collection in arbitrary batches
# from: http://stackoverflow.com/questions/4628290/pairs-from-single-list


def batchwise(coll, size):
    iter_coll = iter(coll)
    return zip(*[iter_coll] * size)

# Returns true if the pair of records are properly pairs


def validpair(read0, read1):
    return read0.name[-1] == "1" and \
        read1.name[-1] == "2" and \
        read0.name[0:-1] == read1.name[0:-1]


# pylint: disable=too-many-locals,too-many-branches
def normalize_by_median(input_filename, htable, args, report_fp=None,
                        report_frequency=100000):

    desired_coverage = args.cutoff
    ksize = htable.ksize()

    # In paired mode we read two records at a time
    batch_size = 1
    if args.paired:
        batch_size = 2

    index = -1
    total = 0
    discarded = 0
    for index, batch in enumerate(batchwise(screed.open(
            input_filename), batch_size)):
        if index > 0 and index % report_frequency == 0:
            print('... kept {kept} of {total} or {perc:2}%'.format(
                kept=total - discarded, total=total,
                perc=int(100. - discarded / float(total) * 100.)))
            print('... in file', input_filename)

            if report_fp:
                print(total, total - discarded, \
                    1. - (discarded / float(total)), file=report_fp)
                report_fp.flush()

        total += batch_size

        # If in paired mode, check that the reads are properly interleaved
        if args.paired:
            if not validpair(batch[0], batch[1]):
                raise IOError('Error: Improperly interleaved pairs \
                    {b0} {b1}'.format(b0=batch[0].name, b1=batch[1].name))

        # Emit the batch of reads if any read passes the filter
        # and all reads are longer than K
        passed_filter = False
        passed_length = True
        for record in batch:
            if len(record.sequence) < ksize:
                passed_length = False
                continue

            seq = record.sequence.replace('N', 'A')
            med, _, _ = htable.get_median_count(seq)

            if med < desired_coverage:
                passed_filter = True
            htable.consume(seq) # consume sequence, either way!

        # Emit records if any passed
        if not (passed_length and passed_filter):
            discarded += batch_size

    return total, discarded


def handle_error(error, input_name):
    print('** ERROR:', error, file=sys.stderr)
    print('** Failed on {name}: '.format(name=input_name), file=sys.stderr)

def get_parser():
    epilog = ("""
    Discard sequences based on whether or not their median k-mer abundance lies
    above a specified cutoff. Kept sequences will be placed in <fileN>.keep.

    Paired end reads will be considered together if :option:`-p` is set. If
    either read will be kept, then both will be kept. This should result in
    keeping (or discarding) each sequencing fragment. This helps with retention
    of repeats, especially.

    With :option:`-s`/:option:`--savegraph`, the k-mer countgraph
    will be saved to the specified file after all sequences have been
    processed. With :option:`-d`, the k-mer countgraph will be
    saved every d files for multifile runs; if :option:`-s` is set,
    the specified name will be used, and if not, the name `backup.ct`
    will be used.  :option:`-l`/:option:`--loadgraph` will load the
    specified k-mer countgraph before processing the specified
    files.  Note that these tables are are in the same format as those
    produced by :program:`load-into-counting.py` and consumed by
    :program:`abundance-dist.py`.

    :option:`-f`/:option:`--fault-tolerant` will force the program to continue
    upon encountering a formatting error in a sequence file; the k-mer counting
    table up to that point will be dumped, and processing will continue on the
    next file.

    Example::

        saturate-by-median.py -k 17 tests/test-data/test-abund-read-2.fa

    Example::

""" "        saturate-by-median.py -p -k 17 tests/test-data/test-abund-read-paired.fa"  # noqa
    """

    Example::

""" "        saturate-by-median.py -k 17 -f tests/test-data/test-error-reads.fq tests/test-data/test-fastq-reads.fq"  # noqa
    """

    Example::

""" "        saturate-by-median.py -k 17 -d 2 -s test.ct tests/test-data/test-abund-read-2.fa tests/test-data/test-fastq-reads")   # noqa
    parser = build_counting_args(
        descr="Do digital normalization (remove mostly redundant sequences)",
        epilog=textwrap.dedent(epilog))
    parser.add_argument('-C', '--cutoff', type=int,
                        default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('-p', '--paired', action='store_true')
    parser.add_argument('-s', '--savegraph', metavar="filename", default='')
    parser.add_argument('-R', '--report',
                        metavar='filename', type=argparse.FileType('w'))
    parser.add_argument('--report-frequency',
                        metavar='report_frequency', default=100000, type=int)
    parser.add_argument('-f', '--fault-tolerant', dest='force',
                        help='continue on next file if read errors are \
                         encountered', action='store_true')
    parser.add_argument('-o', '--out', metavar="filename",
                        dest='single_output_filename',
                        default='', help='only output a single'
                        ' file with the specified filename')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')
    add_loadgraph_args(parser)
    return parser


def main():  # pylint: disable=too-many-branches,too-many-statements
    info('saturate-by-median.py', ['diginorm'])
    parser = sanitize_help(get_parser())
    args = parser.parse_args()

    report_on_config(args)

    report_fp = args.report
    report_frequency = args.report_frequency

    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, False)
    if args.savegraph:
        check_space_for_graph(args, 'countgraph', False)

    # list to save error files along with throwing exceptions
    if args.force:
        corrupt_files = []

    if args.loadgraph:
        print('loading k-mer countgraph from', args.loadgraph)
        htable = khmer.load_countgraph(args.loadgraph)
    else:
        print('making countgraph')
        htable = create_countgraph(args)

    total = 0
    discarded = 0

    for index, input_filename in enumerate(args.input_filenames):
        total_acc = 0
        discarded_acc = 0

        try:
            total_acc, discarded_acc = normalize_by_median(input_filename,
                                                           htable, args,
                                                           report_fp,
                                                           report_frequency)
        except IOError as err:
            handle_error(err, input_filename)
            if not args.force:
                print("NOTE: This can be overridden using the --force"
                      " argument", file=sys.stderr)
                print('** Exiting!', file=sys.stderr)
                sys.exit(1)
            else:
                print('*** Skipping error file, moving on...', file=sys.stderr)
                corrupt_files.append(input_filename)
        else:
            if total_acc == 0 and discarded_acc == 0:
                print('SKIPPED empty file', input_filename)
            else:
                total += total_acc
                discarded += discarded_acc
                print('DONE with {inp}; kept {kept} of {total} or {perc:2}%'\
                    .format(inp=input_filename,
                            kept=total - discarded, total=total,
                            perc=int(100. - discarded / float(total) * 100.)))

    if args.savegraph:
        print('Saving k-mer countgraph through', input_filename)
        print('...saving to', args.savegraph)
        htable.save(args.savegraph)

    # re: threshold, see Zhang et al.,
    # http://arxiv.org/abs/1309.2975
    fp_rate = khmer.calc_expected_collisions(htable, args.force, max_false_pos=.8)
    print('fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate))

    if args.force and len(corrupt_files) > 0:
        print("** WARNING: Finished with errors!", file=sys.stderr)
        print("** I/O Errors occurred in the following files:", file=sys.stderr)
        print("\t", " ".join(corrupt_files), file=sys.stderr)


if __name__ == '__main__':
    main()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
