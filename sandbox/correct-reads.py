#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Semi-streaming error correction.

Output sequences will be placed in 'infile.corr'.

% python scripts/correct-reads.py [ <data1> [ <data2> [ ... ] ] ]

Use -h for parameter help.

TODO: add to sandbox/README.
"""
import sys
import screed
import os
import khmer
import tempfile
import shutil
import textwrap
import argparse

from khmer.khmer_args import (build_counting_args, info, add_loadgraph_args,
                              report_on_config)
from khmer.utils import write_record, write_record_pair, broken_paired_reader
from khmer.kfile import (check_space, check_space_for_graph,
                         check_valid_file_exists)

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2


def correct_sequence(aligner, sequence):
    # align to graph.
    score, graph_alignment, read_alignment, truncated = \
           aligner.align(sequence)

    # next, decide whether or to keep it.
    output_corrected = False
    if not truncated:
        graph_seq = graph_alignment.replace("-", "")
        return True, graph_seq

    return False, sequence


def fix_quality(record):
    if len(record.sequence) < len(record.quality):
        record.quality = record.quality[:len(record.sequence)]

    while len(record.sequence) > len(record.quality):
        record.quality += 'I' # @CTB hack


def get_parser():
    epilog = """
    The output is one file for each input file, <input file>.corr, placed
    in the current directory.  This output contains the input sequences,
    corrected at low-abundance k-mers.

    Note that the output reads will not necessarily be in the same
    order as the reads in the input files. However, read pairs will be
    kept together, in "broken-paired" format; you can use
    ``extract-paired-reads.py`` to extract read pairs and orphans.

    Example::

        correct-reads.py -x 5e7 -k 20 -C 2 data/100k-filtered.fa
    """

    parser = build_counting_args(
        descr='Correct reads using a semi-streaming algorithm.',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('input_filenames', nargs='+')

    parser.add_argument('--cutoff', '-C', type=int,
                        help='k-mers below this abundance are not trusted',
                        default=DEFAULT_CUTOFF)

    parser.add_argument('--normalize-to', '-Z', type=int,
                        help='base cutoff on this median k-mer abundance',
                        default=DEFAULT_NORMALIZE_LIMIT)

    parser.add_argument('-o', '--out', metavar="filename",
                        type=argparse.FileType('w'),
                        default=None, help='only output a single file with '
                        'the specified filename; use a single dash "-" to '
                        'specify that output should go to STDOUT (the '
                        'terminal)')

    parser.add_argument('--variable-coverage', '-V', action='store_true',
                        default=False,
                        help='Only correct sequences that have high coverage.')

    add_loadgraph_args(parser)
    parser.add_argument('-s', '--savegraph', metavar="filename", default='',
                        help='save the k-mer countgraph to disk after all'
                        'reads are loaded.')

    # expert options
    parser.add_argument('--force', default=False, action='store_true')
    parser.add_argument('--ignore-pairs', default=False, action='store_true')
    parser.add_argument('--tempdir', '-T', type=str, default='./')
    parser.add_argument("--theta", dest="bits_theta", type=float, default=1.0)

    return parser


def main():
    info('correct-reads.py', ['streaming'])
    parser = get_parser()
    args = parser.parse_args()

    ###

    if len(set(args.input_filenames)) != len(args.input_filenames):
        print >>sys.stderr, \
            "Error: Cannot input the same filename multiple times."
        sys.exit(1)

    ###

    report_on_config(args)
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savegraph:
        check_space_for_graph(
            args.n_tables * args.min_tablesize, args.force)

    K = args.ksize

    CUTOFF = args.cutoff
    NORMALIZE_LIMIT = args.normalize_to

    if args.loadgraph:
        print >>sys.stderr, 'loading k-mer countgraph from', args.loadgraph
        ct = khmer.load_countinggraph(args.loadgraph)
    else:
        print >>sys.stderr, 'making k-mer countgraph'
        ct = khmer.new_countinggraph(K, args.min_tablesize, args.n_tables)

    tempdir = tempfile.mkdtemp('khmer', 'tmp', args.tempdir)
    print >>sys.stderr, 'created temporary directory %s; ' \
                        'use -T to change location' % tempdir

    aligner = khmer.ReadAligner(ct, args.cutoff, args.bits_theta)

    # ### FIRST PASS ###

    save_pass2_total = 0

    n_bp = 0
    n_reads = 0
    written_bp = 0
    written_reads = 0
    corrected_reads = 0

    pass2list = []
    for filename in args.input_filenames:
        pass2filename = os.path.basename(filename) + '.pass2'
        pass2filename = os.path.join(tempdir, pass2filename)
        if args.out is None:
            corrfp = open(os.path.basename(filename) + '.corr', 'w')
        else:
            corrfp = args.out

        pass2list.append((filename, pass2filename, corrfp))

        screed_iter = screed.open(filename, parse_description=False)
        pass2fp = open(pass2filename, 'w')

        save_pass2 = 0
        n = 0

        paired_iter = broken_paired_reader(screed_iter, min_length=K,
                                           force_single=args.ignore_pairs)
        for n, is_pair, read1, read2 in paired_iter:
            if n % 10000 == 0:
                print >>sys.stderr, '...', n, filename, save_pass2, \
                    n_reads, n_bp, written_reads, written_bp

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
                    is_aligned, new_seq1 = correct_sequence(aligner, seq1)
                    if is_aligned:
                        if new_seq1 != read1.sequence:
                            corrected_reads += 1
                        read1.sequence = new_seq1
                        if hasattr(read1, 'quality'):
                            fix_quality(read1)

                    is_aligned, new_seq2 = correct_sequence(aligner, seq2)
                    if is_aligned:
                        if new_seq2 != read2.sequence:
                            corrected_reads += 1
                        read2.sequence = new_seq2
                        if hasattr(read2, 'quality'):
                            fix_quality(read2)

                    write_record_pair(read1, read2, corrfp)
                    written_reads += 2
                    written_bp += len(read1)
                    written_bp += len(read2)
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
                    is_aligned, new_seq = correct_sequence(aligner, seq)
                    if is_aligned:
                        if new_seq != read1.sequence:
                            corrected_reads += 1
                        read1.sequence = new_seq
                        if hasattr(read1, 'quality'):
                            fix_quality(read1)

                        write_record(read1, corrfp)

                        written_reads += 1
                        written_bp += len(new_seq)

        pass2fp.close()

        print >>sys.stderr, '%s: kept aside %d of %d from first pass, in %s' \
            % (filename, save_pass2, n, filename)
        save_pass2_total += save_pass2

    # ### SECOND PASS. ###

    skipped_n = 0
    skipped_bp = 0
    for _, pass2filename, corrfp in pass2list:
        print >>sys.stderr, ('second pass: looking at sequences kept aside '
                             'in %s') % pass2filename

        # note that for this second pass, we don't care about paired
        # reads - they will be output in the same order they're read in,
        # so pairs will stay together if not orphaned.  This is in contrast
        # to the first loop.

        for n, read in enumerate(screed.open(pass2filename,
                                             parse_description=False)):
            if n % 10000 == 0:
                print >>sys.stderr, '... x 2', n, pass2filename, \
                    written_reads, written_bp

            seq = read.sequence.replace('N', 'A')
            med, _, _ = ct.get_median_count(seq)

            # do we retain low-abundance components unchanged?
            if med < NORMALIZE_LIMIT and args.variable_coverage:
                write_record(read, corrfp)

                written_reads += 1
                written_bp += len(read.sequence)
                skipped_n += 1
                skipped_bp += len(read.sequence)

            # otherwise, examine/correct.
            else:    # med >= NORMALIZE LIMIT or not args.variable_coverage
                is_aligned, new_seq = correct_sequence(aligner, seq)
                if is_aligned:
                    if new_seq != read.sequence:
                        corrected_reads += 1
                    read.sequence = new_seq
                    if hasattr(read, 'quality'):
                        fix_quality(read)
                    write_record(read, corrfp)

                    written_reads += 1
                    written_bp += len(new_seq)

        print >>sys.stderr, 'removing %s' % pass2filename
        os.unlink(pass2filename)

    print >>sys.stderr, 'removing temp directory & contents (%s)' % tempdir
    shutil.rmtree(tempdir)

    n_passes = 1.0 + (float(save_pass2_total) / n_reads)
    percent_reads_corrected = float(corrected_reads +
                                    (n_reads - written_reads)) /\
        n_reads * 100.0

    print >>sys.stderr, 'read %d reads, %d bp' % (n_reads, n_bp,)
    print >>sys.stderr, 'wrote %d reads, %d bp' % (written_reads, written_bp,)
    print >>sys.stderr, 'looked at %d reads twice (%.2f passes)' % \
        (save_pass2_total, n_passes)
    print >>sys.stderr, 'removed %d reads and corrected %d reads (%.2f%%)' % \
        (n_reads - written_reads, corrected_reads, percent_reads_corrected)
    print >>sys.stderr, 'removed %.2f%% of bases (%d total)' % \
        ((1 - (written_bp / float(n_bp))) * 100.0, n_bp - written_bp)

    if args.variable_coverage:
        percent_reads_hicov = 100.0 * float(n_reads - skipped_n) / n_reads
        print >>sys.stderr, '%d reads were high coverage (%.2f%%);' % \
            (n_reads - skipped_n, percent_reads_hicov)
        print >>sys.stderr, ('skipped %d reads/%d bases because of low'
                             'coverage') % (skipped_n, skipped_bp)

    fp_rate = \
        khmer.calc_expected_collisions(ct, args.force, max_false_pos=.8)
    # for max_false_pos see Zhang et al., http://arxiv.org/abs/1309.2975
    print >>sys.stderr, \
        'fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate)

    print >>sys.stderr, 'output in *.corr'

    if args.savegraph:
        print >>sys.stderr, "Saving k-mer countgraph to", args.savegraph
        ct.save(args.savegraph)


if __name__ == '__main__':
    main()
