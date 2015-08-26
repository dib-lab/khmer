#! /usr/bin/env python
# Copyright (C) 2012-2015, Michigan State University.
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
# pylint: disable=invalid-name,missing-docstring,no-member
"""
Trim sequences at k-mers of the given abundance, using a streaming algorithm.

Output sequences will be placed in 'infile.abundtrim'.

% python scripts/trim-low-abund.py [ <data1> [ <data2> [ ... ] ] ]

Use -h for parameter help.
"""
from __future__ import print_function
import sys
import screed
import os
import khmer
import tempfile
import shutil
import textwrap
import argparse

from screed import Record
from khmer import khmer_args

from khmer.khmer_args import (build_counting_args, info, add_loadgraph_args,
                              report_on_config, calculate_graphsize,
                              sanitize_help)
from khmer.utils import write_record, write_record_pair, broken_paired_reader
from khmer.kfile import (check_space, check_space_for_graph,
                         check_valid_file_exists, add_output_compression_type,
                         get_file_writer)

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2


def trim_record(read, trim_at):
    new_read = Record()
    new_read.name = read.name
    new_read.sequence = read.sequence[:trim_at]
    if hasattr(read, 'quality'):
        new_read.quality = read.quality[:trim_at]

    return new_read


def get_parser():
    epilog = """\
    The output is one file for each input file, ``<input file>.abundtrim``,
    placed in the current directory.  This output contains the input sequences
    trimmed at low-abundance k-mers.

    The :option:`-V`/:option:`--variable-coverage` parameter will, if
    specified, prevent elimination of low-abundance reads by only trimming
    low-abundance k-mers from high-abundance reads; use this for
    non-genomic data sets that may have variable coverage.

    Note that the output reads will not necessarily be in the same order
    as the reads in the input files; if this is an important consideration,
    use :program:`load-into-counting.py` and :program:`filter-abund.py`.
    However, read pairs will be kept together, in "broken-paired" format; you
    can use :program:`extract-paired-reads.py` to extract read pairs and
    orphans.

    Example::

        trim-low-abund.py -x 5e7 -k 20 -C 2 data/100k-filtered.fa
    """

    parser = build_counting_args(
        descr='Trim low-abundance k-mers using a streaming algorithm.',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('input_filenames', nargs='+')

    parser.add_argument('--cutoff', '-C', type=int,
                        help='remove k-mers below this abundance',
                        default=DEFAULT_CUTOFF)

    parser.add_argument('--normalize-to', '-Z', type=int,
                        help='base cutoff on this median k-mer abundance',
                        default=DEFAULT_NORMALIZE_LIMIT)

    parser.add_argument('-o', '--output', metavar="output_filename",
                        type=argparse.FileType('wb'),
                        help='only output a single file with '
                        'the specified filename; use a single dash "-" to '
                        'specify that output should go to STDOUT (the '
                        'terminal)')

    parser.add_argument('--variable-coverage', '-V', action='store_true',
                        default=False,
                        help='Only trim low-abundance k-mers from sequences '
                        'that have high coverage.')

    add_loadgraph_args(parser)
    parser.add_argument('-s', '--savegraph', metavar="filename", default='',
                        help='save the k-mer countgraph to disk after all'
                        'reads are loaded.')

    # expert options
    parser.add_argument('--force', default=False, action='store_true')
    parser.add_argument('--ignore-pairs', default=False, action='store_true')
    parser.add_argument('--tempdir', '-T', type=str, default='./')
    add_output_compression_type(parser)

    return parser


def main():
    info('trim-low-abund.py', ['streaming'])
    parser = sanitize_help(get_parser())
    args = parser.parse_args()

    ###

    if len(set(args.input_filenames)) != len(args.input_filenames):
        print("Error: Cannot input the same filename multiple times.",
              file=sys.stderr)
        sys.exit(1)

    ###

    report_on_config(args)
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savegraph:
        graphsize = calculate_graphsize(args, 'countgraph')
        check_space_for_graph(args.savegraph, graphsize, args.force)

    if ('-' in args.input_filenames or '/dev/stdin' in args.input_filenames) \
       and not args.output:
        print("Accepting input from stdin; output filename must "
              "be provided with -o.", file=sys.stderr)
        sys.exit(1)

    if args.loadgraph:
        print('loading countgraph from', args.loadgraph, file=sys.stderr)
        ct = khmer.load_countgraph(args.loadgraph)
    else:
        print('making countgraph', file=sys.stderr)
        ct = khmer_args.create_countgraph(args)

    K = ct.ksize()
    CUTOFF = args.cutoff
    NORMALIZE_LIMIT = args.normalize_to

    tempdir = tempfile.mkdtemp('khmer', 'tmp', args.tempdir)
    print('created temporary directory %s; '
          'use -T to change location' % tempdir, file=sys.stderr)

    # ### FIRST PASS ###

    save_pass2_total = 0

    n_bp = 0
    n_reads = 0
    written_bp = 0
    written_reads = 0
    trimmed_reads = 0

    pass2list = []
    for filename in args.input_filenames:
        pass2filename = os.path.basename(filename) + '.pass2'
        pass2filename = os.path.join(tempdir, pass2filename)
        if args.output is None:
            trimfp = get_file_writer(open(os.path.basename(filename) +
                                          '.abundtrim', 'wb'),
                                     args.gzip, args.bzip)
        else:
            trimfp = get_file_writer(args.output, args.gzip, args.bzip)

        pass2list.append((filename, pass2filename, trimfp))

        screed_iter = screed.open(filename)
        pass2fp = open(pass2filename, 'w')

        save_pass2 = 0
        n = 0

        paired_iter = broken_paired_reader(screed_iter, min_length=K,
                                           force_single=args.ignore_pairs)
        for n, is_pair, read1, read2 in paired_iter:
            if n % 10000 == 0:
                print('...', n, filename, save_pass2, n_reads, n_bp,
                      written_reads, written_bp, file=sys.stderr)

            # we want to track paired reads here, to make sure that pairs
            # are not split between first pass and second pass.

            if is_pair:
                n_reads += 2
                n_bp += len(read1.sequence) + len(read2.sequence)

                seq1 = read1.sequence.replace('N', 'A')
                seq2 = read2.sequence.replace('N', 'A')

                med1, _, _ = ct.get_median_count(seq1)
                med2, _, _ = ct.get_median_count(seq2)

                if med1 < NORMALIZE_LIMIT or med2 < NORMALIZE_LIMIT:
                    ct.consume(seq1)
                    ct.consume(seq2)
                    write_record_pair(read1, read2, pass2fp)
                    save_pass2 += 2
                else:
                    _, trim_at1 = ct.trim_on_abundance(seq1, CUTOFF)
                    _, trim_at2 = ct.trim_on_abundance(seq2, CUTOFF)

                    if trim_at1 >= K:
                        read1 = trim_record(read1, trim_at1)

                    if trim_at2 >= K:
                        read2 = trim_record(read2, trim_at2)

                    if trim_at1 != len(seq1):
                        trimmed_reads += 1
                    if trim_at2 != len(seq2):
                        trimmed_reads += 1

                    write_record_pair(read1, read2, trimfp)
                    written_reads += 2
                    written_bp += trim_at1 + trim_at2
            else:
                n_reads += 1
                n_bp += len(read1.sequence)

                seq = read1.sequence.replace('N', 'A')

                med, _, _ = ct.get_median_count(seq)

                # has this portion of the graph saturated? if not,
                # consume & save => pass2.
                if med < NORMALIZE_LIMIT:
                    ct.consume(seq)
                    write_record(read1, pass2fp)
                    save_pass2 += 1
                else:                       # trim!!
                    _, trim_at = ct.trim_on_abundance(seq, CUTOFF)
                    if trim_at >= K:
                        new_read = trim_record(read1, trim_at)
                        write_record(new_read, trimfp)

                        written_reads += 1
                        written_bp += trim_at

                        if trim_at != len(read1.sequence):
                            trimmed_reads += 1

        pass2fp.close()

        print('%s: kept aside %d of %d from first pass, in %s' %
              (filename, save_pass2, n, filename),
              file=sys.stderr)
        save_pass2_total += save_pass2

    # ### SECOND PASS. ###

    skipped_n = 0
    skipped_bp = 0
    for _, pass2filename, trimfp in pass2list:
        print('second pass: looking at sequences kept aside in %s' %
              pass2filename,
              file=sys.stderr)

        # note that for this second pass, we don't care about paired
        # reads - they will be output in the same order they're read in,
        # so pairs will stay together if not orphaned.  This is in contrast
        # to the first loop.

        for n, read in enumerate(screed.open(pass2filename)):
            if n % 10000 == 0:
                print('... x 2', n, pass2filename,
                      written_reads, written_bp, file=sys.stderr)

            seq = read.sequence.replace('N', 'A')
            med, _, _ = ct.get_median_count(seq)

            # do we retain low-abundance components unchanged?
            if med < NORMALIZE_LIMIT and args.variable_coverage:
                write_record(read, trimfp)

                written_reads += 1
                written_bp += len(read.sequence)
                skipped_n += 1
                skipped_bp += len(read.sequence)

            # otherwise, examine/trim/truncate.
            else:    # med >= NORMALIZE LIMIT or not args.variable_coverage
                _, trim_at = ct.trim_on_abundance(seq, CUTOFF)
                if trim_at >= K:
                    new_read = trim_record(read, trim_at)
                    write_record(new_read, trimfp)

                    written_reads += 1
                    written_bp += trim_at

                    if trim_at != len(read.sequence):
                        trimmed_reads += 1

        print('removing %s' % pass2filename, file=sys.stderr)
        os.unlink(pass2filename)

    print('removing temp directory & contents (%s)' % tempdir, file=sys.stderr)
    shutil.rmtree(tempdir)

    n_passes = 1.0 + (float(save_pass2_total) / n_reads)
    percent_reads_trimmed = float(trimmed_reads + (n_reads - written_reads)) /\
        n_reads * 100.0

    print('read %d reads, %d bp' % (n_reads, n_bp,), file=sys.stderr)
    print('wrote %d reads, %d bp' % (written_reads, written_bp,),
          file=sys.stderr)
    print('looked at %d reads twice (%.2f passes)' % (save_pass2_total,
                                                      n_passes),
          file=sys.stderr)
    print('removed %d reads and trimmed %d reads (%.2f%%)' %
          (n_reads - written_reads, trimmed_reads, percent_reads_trimmed),
          file=sys.stderr)
    print('trimmed or removed %.2f%% of bases (%d total)' %
          ((1 - (written_bp / float(n_bp))) * 100.0, n_bp - written_bp),
          file=sys.stderr)

    if args.variable_coverage:
        percent_reads_hicov = 100.0 * float(n_reads - skipped_n) / n_reads
        print('%d reads were high coverage (%.2f%%);' % (n_reads - skipped_n,
                                                         percent_reads_hicov),
              file=sys.stderr)
        print('skipped %d reads/%d bases because of low coverage' %
              (skipped_n, skipped_bp),
              file=sys.stderr)

    fp_rate = \
        khmer.calc_expected_collisions(ct, args.force, max_false_pos=.8)
    # for max_false_pos see Zhang et al., http://arxiv.org/abs/1309.2975
    print('fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate),
          file=sys.stderr)

    print('output in *.abundtrim', file=sys.stderr)

    if args.savegraph:
        print("Saving k-mer countgraph to",
              args.savegraph, file=sys.stderr)
        ct.save(args.savegraph)


if __name__ == '__main__':
    main()
