#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
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
                              sanitize_epilog)
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
    epilog = """
    The output is one file for each input file, <input file>.abundtrim, placed
    in the current directory.  This output contains the input sequences
    trimmed at low-abundance k-mers.

    The ``-V/--variable-coverage`` parameter will, if specified,
    prevent elimination of low-abundance reads by only trimming
    low-abundance k-mers from high-abundance reads; use this for
    non-genomic data sets that may have variable coverage.

    Note that the output reads will not necessarily be in the same order
    as the reads in the input files; if this is an important consideration,
    use ``load-into-countgraph.py`` and ``filter-abund.py``.  However, read
    pairs will be kept together, in "broken-paired" format; you can use
    ``extract-paired-reads.py`` to extract read pairs and orphans.

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

    parser.add_argument('--diginorm', default=False, action='store_true')
    parser.add_argument('--single-pass', default=False, action='store_true')

    return parser


class Trimmer(object):
    def __init__(self, graph, do_trim_low_abund, cutoff, normalize_limit):
        self.graph = graph
        self.do_trim_low_abund = do_trim_low_abund
        self.cutoff = cutoff
        self.normalize_limit = normalize_limit

        self.n_reads = 0
        self.n_bp = 0
        self.trimmed_reads = 0
        self.n_saved = 0
        self.n_skipped = 0
        self.bp_skipped = 0

    def __call__(self, reader, saver):
        graph = self.graph
        NORMALIZE_LIMIT = self.normalize_limit
        CUTOFF = self.cutoff
        K = graph.ksize()

        for n, is_pair, read1, read2 in reader:
            examine = []

            if is_pair:
                reads = (read1, read2)
            else:
                reads = (read1,)

            # clean up the sequences for examination.
            for read in reads:
                r = read.sequence.replace('N', 'A')
                examine.append(r)

                self.n_reads += 1
                self.n_bp += len(r)

            # find out if they are estimated to have low coverage
            is_low_coverage = False
            for r in examine:
                med, _, _ = graph.get_median_count(r)
                if med < NORMALIZE_LIMIT:
                    is_low_coverage = True
                    break

            # if either of the sequences are low coverage & we have a 'saver',
            # keep both for 2nd pass
            if is_low_coverage and saver:
                for read, seq in zip(reads, examine):
                    graph.consume(seq)
                    write_record(read, saver)
                    self.n_saved += 1

            # if they're low coverage, and we don't want to trim low coverage,
            # write them out as is.
            elif is_low_coverage and not self.do_trim_low_abund:
                for read in reads:
                    self.n_skipped += 1
                    self.bp_skipped += len(read)
                    yield read
            # otherwise, trim them if they should be trimmed, THEN write 'em
            else:
                assert (not is_low_coverage or self.do_trim_low_abund)
                for read, seq in zip(reads, examine):
                    _, trim_at = graph.trim_on_abundance(seq, CUTOFF)

                    # too short after trimming? eliminate read.
                    if trim_at < K:
                        continue

                    # will trim? do so.
                    if trim_at != len(seq):
                        self.trimmed_reads += 1
                        read = trim_record(read, trim_at)

                    # save
                    yield read


def main():
    info('trim-low-abund.py', ['streaming'])
    parser = sanitize_epilog(get_parser())
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

    written_bp = 0
    written_reads = 0

    trimmer = Trimmer(ct, not args.variable_coverage, args.cutoff,
                      args.normalize_to)

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

        paired_iter = broken_paired_reader(screed_iter, min_length=K,
                                           force_single=args.ignore_pairs)

        n_start = trimmer.n_reads
        save_start = trimmer.n_saved
        for read in trimmer(paired_iter, pass2fp):
            if (trimmer.n_reads - n_start) % 10000 == 0:
                print('...', n, filename, trimmer.n_saved,
                      trimmer.n_reads, trimmer.n_bp,
                      written_reads, written_bp, file=sys.stderr)

            write_record(read, trimfp)
            written_bp += len(read)
            written_reads += 1
        pass2fp.close()

        print('%s: kept aside %d of %d from first pass, in %s' %
              (filename,
               trimmer.n_saved - save_start, trimmer.n_reads - n_start,
               filename),
              file=sys.stderr)

    # first pass goes across all the data, so record relevant stats...
    n_reads = trimmer.n_reads           
    n_bp = trimmer.n_bp
    n_skipped = trimmer.n_skipped
    bp_skipped = trimmer.bp_skipped
    save_pass2_total = trimmer.n_saved

    # ### SECOND PASS. ###

    # nothing should have been skipped yet!
    assert trimmer.n_skipped == 0
    assert trimmer.bp_skipped == 0

    if args.single_pass:
        pass2list = []

    # go back through all the files again.
    for _, pass2filename, trimfp in pass2list:
        print('second pass: looking at sequences kept aside in %s' %
              pass2filename,
              file=sys.stderr)

        # note that for this second pass, we don't care about paired
        # reads - they will be output in the same order they're read in,
        # so pairs will stay together if not orphaned.  This is in contrast
        # to the first loop.  Hence, force_single=True below.

        screed_iter = screed.open(pass2filename, parse_description=False)
        paired_iter = broken_paired_reader(screed_iter, min_length=K,
                                           force_single=True)

        for read in trimmer(paired_iter, None):
            if (trimmer.n_reads - n_start) % 10000 == 0:
                print('... x 2', trimmer.n_reads - n_start,
                      pass2filename, trimmer.n_saved,
                      trimmer.n_reads, trimmer.n_bp,
                      written_reads, written_bp, file=sys.stderr)

            write_record(read, trimfp)
            written_reads += 1
            written_bp += len(read)

        print('removing %s' % pass2filename, file=sys.stderr)
        os.unlink(pass2filename)

    print('removing temp directory & contents (%s)' % tempdir, file=sys.stderr)
    shutil.rmtree(tempdir)

    trimmed_reads = trimmer.trimmed_reads

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
        percent_reads_hicov = 100.0 * float(n_reads - n_skipped) / n_reads
        print('%d reads were high coverage (%.2f%%);' % (n_reads - n_skipped,
                                                         percent_reads_hicov),
              file=sys.stderr)
        print('skipped %d reads/%d bases because of low coverage' %
              (n_skipped, bp_skipped),
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
