#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring
"""
Eliminate surplus reads.

Eliminate reads with median k-mer abundance higher than
DESIRED_COVERAGE.  Output sequences will be placed in 'infile.keep', with the
option to output to STDOUT.

% python scripts/normalize-by-median.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import screed
import os
import khmer
import textwrap
from khmer import khmer_args
from contextlib import contextmanager
from oxli import functions as oxutils
from khmer.khmer_args import (build_counting_args, add_loadhash_args,
                              report_on_config, info, calculate_tablesize)
import argparse
from khmer.kfile import (check_space, check_space_for_hashtable,
                         check_valid_file_exists)
from khmer.utils import write_record, broken_paired_reader
from khmer.khmer_logger import (configure_logging, log_info, log_error)


DEFAULT_DESIRED_COVERAGE = 20


class WithDiagnostics(object):
    """
    Generator/context manager to do boilerplate output of statistics.

    uses a Normalizer object.
    """

    def __init__(self, norm, report_fp=None, report_frequency=100000):
        self.norm = norm
        self.report_fp = report_fp
        if report_fp:
            report_fp.write('total,kept,f_kept\n')

        self.total = 0
        self.kept = 0

        self.report_frequency = report_frequency
        self.next_report_at = self.report_frequency
        self.last_report_at = self.report_frequency

    def __call__(self, reader, ifilename):
        norm = self.norm
        report_fp = self.report_fp

        reads_start = self.total
        total = self.total
        kept = self.kept

        try:
            for _, is_paired, read0, read1 in reader:
                if is_paired:
                    total += 2
                else:
                    total += 1

                # do diginorm
                for record in norm(is_paired, read0, read1):
                    kept += 1
                    yield record

                # report!
                if total >= self.next_report_at:
                    self.next_report_at += self.report_frequency
                    self.last_report_at = total

                    perc_kept = kept / float(total)

                    log_info('... kept {kept} of {tot} or {perc_kept:.1%} so'
                             'far', kept=kept, tot=total,
                             perc_kept=perc_kept)
                    log_info('... in file ' + ifilename)

                    if report_fp:
                        print("{total},{kept},{f_kept:.4}"
                              .format(total=total, f_kept=perc_kept,
                                      kept=kept),
                              file=report_fp)
                        report_fp.flush()
        finally:
            self.total = total
            self.kept = kept

        # per file diagnostic output
        if total == reads_start:
            log_info('SKIPPED empty file ' + ifilename)
        else:
            perc_kept = kept / float(total)

            log_info('DONE with {inp}; kept {kept} of {total} or '
                     '{perc_kept:.1%}', inp=ifilename, kept=kept, total=total,
                     perc_kept=perc_kept)

        # make sure there's at least one report per file, at the end of each
        # file.
        if report_fp and total != self.last_report_at:
            perc_kept = kept / float(total)

            print("{total},{kept},{f_kept:.4}"
                  .format(total=total, f_kept=perc_kept, kept=kept),
                  file=report_fp)
            report_fp.flush()


class Normalizer(object):

    """Digital normalization algorithm."""

    def __init__(self, desired_coverage, htable):
        self.htable = htable
        self.desired_coverage = desired_coverage

    def __call__(self, is_paired, read0, read1):
        """
        Actually does digital normalization - the core algorithm.

        * get one (unpaired) or two (paired) reads;
        * sanitize the sequences (convert Ns to As);
        * get the median k-mer count of one/both reads;
        * if any read's median k-mer count is below desired coverage, keep all;
        * consume and yield kept reads.
        """
        desired_coverage = self.desired_coverage

        passed_filter = False

        batch = []
        batch.append(read0)
        if read1 is not None:
            batch.append(read1)

        for record in batch:
            seq = record.sequence.replace('N', 'A')
            if not self.htable.median_at_least(seq, desired_coverage):
                passed_filter = True

        if passed_filter:
            for record in batch:
                seq = record.sequence.replace('N', 'A')
                self.htable.consume(seq)
                yield record


@contextmanager
def catch_io_errors(ifile, out, single_out, force, corrupt_files):
    """Context manager to do boilerplate handling of IOErrors."""
    try:
        yield
    except (IOError, OSError, ValueError) as error:
        log_error('** ERROR: ' + str(error))
        log_error('** Failed on {name}: ', name=ifile)
        if not single_out:
            os.remove(out.name)
        if not force:
            log_error('** Exiting!')

            sys.exit(1)
        else:
            log_error('*** Skipping error file, moving on...')
            corrupt_files.append(ifile)


def get_parser():
    epilog = ("""
    Discard sequences based on whether or not their median k-mer abundance lies
    above a specified cutoff. Kept sequences will be placed in <fileN>.keep.

    By default, paired end reads will be considered together; if
    either read should be kept, both will be kept. (This keeps both
    reads from a fragment, and helps with retention of repeats.)
    Unpaired reads are treated individually.

    If :option:`-p`/`--paired` is set, then proper pairing is required
    and the script will exit on unpaired reads, although
    :option:`--unpaired-reads` can be used to supply a file of orphan
    reads to be read after the paired reads.

    :option:`--force-single` will ignore all pairing information and treat
    reads individually.

    With :option:`-s`/:option:`--savetable`, the k-mer counting table
    will be saved to the specified file after all sequences have been
    processed. :option:`-l`/:option:`--loadtable` will load the
    specified k-mer counting table before processing the specified
    files.  Note that these tables are are in the same format as those
    produced by :program:`load-into-counting.py` and consumed by
    :program:`abundance-dist.py`.

    To append reads to an output file (rather than overwriting it), send output
    to STDOUT with `--out -` and use UNIX file redirection syntax (`>>`) to
    append to the file.

    Example::

        normalize-by-median.py -k 17 tests/test-data/test-abund-read-2.fa

    Example::

""" "        normalize-by-median.py -p -k 17 tests/test-data/test-abund-read-paired.fa"  # noqa
    """

    Example::

""" "        normalize-by-median.py -p -k 17 -o - tests/test-data/paired.fq >> appended-output.fq"  # noqa
    """

    Example::

""" "        normalize-by-median.py -k 17 -f tests/test-data/test-error-reads.fq tests/test-data/test-fastq-reads.fq"  # noqa
    """

    Example::

""" "        normalize-by-median.py -k 17 -d 2 -s test.ct tests/test-data/test-abund-read-2.fa tests/test-data/test-fastq-reads")   # noqa
    parser = build_counting_args(
        descr="Do digital normalization (remove mostly redundant sequences)",
        epilog=textwrap.dedent(epilog))
    parser.add_argument('-C', '--cutoff', type=int,
                        default=DEFAULT_DESIRED_COVERAGE)
    parser.add_argument('-p', '--paired', action='store_true',
                        help='require that all sequences be properly paired')
    parser.add_argument('--force-single', dest='force_single',
                        action='store_true',
                        help='treat all sequences as single-ended/unpaired')
    parser.add_argument('-u', '--unpaired-reads',
                        metavar="unpaired_reads_filename",
                        help='include a file of unpaired reads to which '
                        '-p/--paired does not apply.')
    parser.add_argument('-s', '--savetable', metavar="filename", default='',
                        help='save the k-mer counting table to disk after all'
                        'reads are loaded.')
    parser.add_argument('-R', '--report',
                        metavar='report_filename', type=argparse.FileType('w'))
    parser.add_argument('--report-frequency',
                        metavar='report_frequency', type=int,
                        default=100000)
    parser.add_argument('-f', '--force', dest='force',
                        help='continue past file reading errors',
                        action='store_true')
    parser.add_argument('-o', '--out', metavar="filename",
                        dest='single_output_file',
                        type=argparse.FileType('w'),
                        default=None, help='only output a single file with '
                        'the specified filename; use a single dash "-" to '
                        'specify that output should go to STDOUT (the '
                        'terminal)')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')
    add_loadhash_args(parser)
    return parser


def main():  # pylint: disable=too-many-branches,too-many-statements

    args = get_parser().parse_args()

    info('normalize-by-median.py', ['diginorm'])
    report_on_config(args)

    configure_logging(args.quiet)

    report_fp = args.report
    force_single = args.force_single

    # if estimates of k-mer numbers are given, choose table sizes.
    args = oxutils.do_sanity_checking(args, 0.1)

    # check for similar filenames
    # if we're using a single output file only check for identical filenames
    # otherwise, check for identical BASE names as well.
    filenames = []
    basenames = []
    for pathfilename in args.input_filenames:
        filenames.append(pathfilename)
        if args.single_output_file:
            continue  # nothing more to worry about

        basename = os.path.basename(pathfilename)
        if basename in basenames:
            log_error('ERROR: Duplicate filename--Cannot handle this!')
            log_error('** Exiting!')
            sys.exit(1)

        basenames.append(basename)

    # check that files exist and there is sufficient output disk space.
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savetable:
        tablesize = calculate_tablesize(args, 'countgraph')
        check_space_for_hashtable(args.savetable, tablesize, args.force)

    # load or create counting table.
    if args.loadtable:
        log_info('loading k-mer counting table from ' + args.loadtable)
        htable = khmer.load_counting_hash(args.loadtable)
        if args.unique_kmers != 0:
            log_info('Warning: You have specified the number of unique kmers'
                     ' but are loading a precreated counting table --'
                     ' counting table size will NOT be set automatically.')
    else:
        log_info('making countgraph')
        htable = khmer_args.create_countgraph(args)

    # create an object to handle diginorm of all files
    norm = Normalizer(args.cutoff, htable)
    with_diagnostics = WithDiagnostics(norm, report_fp, args.report_frequency)

    # make a list of all filenames and if they're paired or not;
    # if we don't know if they're paired, default to allowing but not
    # forcing pairing.
    files = []
    for element in filenames:
        files.append([element, args.paired])
    if args.unpaired_reads:
        files.append([args.unpaired_reads, False])

    corrupt_files = []
    outfp = None
    output_name = None

    if args.single_output_file:
        if args.single_output_file is sys.stdout:
            output_name = '/dev/stdout'
        else:
            output_name = args.single_output_file.name
        outfp = args.single_output_file
    else:
        if '-' in filenames or '/dev/stdin' in filenames:
            print("Accepting input from stdin; output filename must "
                  "be provided with '-o'.", file=sys.stderr)
            sys.exit(1)

    #
    # main loop: iterate over all files given, do diginorm.
    #

    for filename, require_paired in files:
        if not args.single_output_file:
            output_name = os.path.basename(filename) + '.keep'
            outfp = open(output_name, 'w')

        # failsafe context manager in case an input file breaks
        with catch_io_errors(filename, outfp, args.single_output_file,
                             args.force, corrupt_files):

            screed_iter = screed.open(filename)
            reader = broken_paired_reader(screed_iter, min_length=args.ksize,
                                          force_single=force_single,
                                          require_paired=require_paired)

            # actually do diginorm
            for record in with_diagnostics(reader, filename):
                if record is not None:
                    write_record(record, outfp)

            log_info('output in ' + output_name)
            if output_name is not '/dev/stdout':
                outfp.close()

    # finished - print out some diagnostics.

    log_info('Total number of unique k-mers: {umers}',
             umers=htable.n_unique_kmers())

    if args.savetable:
        log_info('...saving to ' + args.savetable)
        htable.save(args.savetable)

    fp_rate = \
        khmer.calc_expected_collisions(htable, False, max_false_pos=.8)
    # for max_false_pos see Zhang et al., http://arxiv.org/abs/1309.2975

    log_info('fp rate estimated to be {fpr:1.3f}', fpr=fp_rate)

    if args.force and len(corrupt_files) > 0:
        log_error("** WARNING: Finished with errors!")
        log_error("** I/O Errors occurred in the following files:")
        log_error("\t" + " ".join(corrupt_files))


if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
