#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Trim sequences at k-mers of the given abundance, using a streaming algorithm.
Output sequences will be placed in 'infile.abundtrim'.

% python scripts/trim-low-abund.py [ <data1> [ <data2> [ ... ] ] ]

Use -h for parameter help.

TODO: load/save counting table.
TODO: reference appropriate preprint.
"""
import sys
import screed
import os
import khmer
import tempfile
import shutil
import textwrap

from screed.screedRecord import _screed_record_dict
from khmer.khmer_args import build_counting_args
from khmer.utils import (write_record, write_record_pair, broken_paired_reader)


DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2

# see Zhang et al., http://arxiv.org/abs/1309.2975
MAX_FALSE_POSITIVE_RATE = 0.8


def trim_record(read, trim_at):
    new_read = _screed_record_dict()
    new_read.name = read.name
    new_read.sequence = read.sequence[:trim_at]
    if hasattr(read, 'accuracy'):
        new_read.accuracy = read.accuracy[:trim_at]

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
    use ``load-into-counting.py`` and ``filter-abund.py``.  However, read
    pairs will be kept together, in "broken-paired" format; you can use
    ``extract-paired-reads.py`` to extract read pairs and orphans.

    Example::

        trim-low-abund.py -x 5e7 -k 20 -C 2 data/100k-filtered.fa
    """

    parser = build_counting_args(
        descr='Trim low-abundance k-mers using a streaming algorithm.',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('--cutoff', '-C', type=int, dest='cutoff',
                        help='remove k-mers below this abundance',
                        default=DEFAULT_CUTOFF)

    parser.add_argument('--normalize-to', '-Z', type=int, dest='normalize_to',
                        help='base cutoff on this median k-mer abundance',
                        default=DEFAULT_NORMALIZE_LIMIT)

    parser.add_argument('--variable-coverage', '-V', action='store_true',
                        dest='variable_coverage', default=False,
                        help='Only trim low-abundance k-mers from sequences '
                        'that have high coverage.')
    parser.add_argument('--tempdir', '-T', type=str, dest='tempdir',
                        default='./')

    parser.add_argument('input_filenames', nargs='+')

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    ###

    if len(set(args.input_filenames)) != len(args.input_filenames):
        print >>sys.stderr, \
            "Error: Cannot input the same filename multiple times."
        sys.exit(1)

    ###

    K = args.ksize

    CUTOFF = args.cutoff
    NORMALIZE_LIMIT = args.normalize_to

    print >>sys.stderr, 'making counting table'
    ct = khmer.new_counting_hash(K, args.min_tablesize, args.n_tables)

    tempdir = tempfile.mkdtemp('khmer', 'tmp', args.tempdir)
    print >>sys.stderr, 'created temporary directory %s; ' \
                        'use -T to change location' % tempdir

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
        trimfilename = os.path.basename(filename) + '.abundtrim'

        pass2list.append((filename, pass2filename, trimfilename))

        screed_iter = screed.open(filename)
        pass2fp = open(pass2filename, 'w')
        trimfp = open(trimfilename, 'w')

        save_pass2 = 0
        for n, is_pair, read1, read2 in broken_paired_reader(screed_iter):
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
        trimfp.close()

        print '%s: kept aside %d of %d from first pass, in %s' % \
              (filename, save_pass2, n, filename)
        save_pass2_total += save_pass2

    # ### SECOND PASS. ###

    skipped_n = 0
    skipped_bp = 0
    for orig_filename, pass2filename, trimfilename in pass2list:
        print 'second pass: looking at sequences kept aside in %s' % \
              pass2filename

        # note that for this second pass, we don't care about paired
        # reads - they will be output in the same order they're read in,
        # so pairs will stay together if not orphaned.  This is in contrast
        # to the first loop.

        trimfp = open(trimfilename, 'a')
        for n, read in enumerate(screed.open(pass2filename)):
            if n % 10000 == 0:
                print >>sys.stderr, '... x 2', n, pass2filename, \
                    written_reads, written_bp

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

        print >>sys.stderr, 'removing %s' % pass2filename
        os.unlink(pass2filename)

    print >>sys.stderr, 'removing temp directory & contents (%s)' % tempdir
    shutil.rmtree(tempdir)

    print 'read %d reads, %d bp' % (n_reads, n_bp,)
    print 'wrote %d reads, %d bp' % (written_reads, written_bp,)
    print 'removed %d reads and trimmed %d reads' % (n_reads - written_reads,
                                                     trimmed_reads,)
    print 'looked at %d reads twice' % (save_pass2_total,)
    print 'trimmed or removed %.2f%% of bases (%d total)' % \
        ((1 - (written_bp / float(n_bp))) * 100., n_bp - written_bp)

    if args.variable_coverage:
        print 'skipped %d reads/%d bases because of low coverage' % \
              (skipped_n, skipped_bp)

    fp_rate = khmer.calc_expected_collisions(ct)
    print >>sys.stderr, \
        'fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate)

    if fp_rate > MAX_FALSE_POSITIVE_RATE:
        print >> sys.stderr, "**"
        print >> sys.stderr, ("** ERROR: the k-mer counting table is too small"
                              " for this data set. Increase tablesize/# "
                              "tables.")
        print >> sys.stderr, "**"
        print >> sys.stderr, "** Do not use these results!!"
        sys.exit(1)

    print 'output in *.abundtrim'


if __name__ == '__main__':
    main()
