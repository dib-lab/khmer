#! /usr/bin/env python2
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
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
from itertools import izip
from contextlib import contextmanager

from khmer.khmer_args import (build_counting_args, add_loadhash_args,
                              report_on_config, info)
import argparse
from khmer.kfile import (check_space, check_space_for_hashtable,
                         check_valid_file_exists)
from khmer.utils import write_record, check_is_pair
DEFAULT_DESIRED_COVERAGE = 10

# Iterate a collection in arbitrary batches
# from: http://stackoverflow.com/questions/4628290/pairs-from-single-list


def batchwise(coll, size):
    iter_coll = iter(coll)
    return izip(*[iter_coll] * size)

# Returns true if the pair of records are properly pairs


def WithDiagnostics(ifile, batch_size, fp, paired, norm):
    """
    Generator/context manager to do boilerplate output of statistics while
    normalizing data. Also checks for properly paired data.
    """

    index = 0

    for index, batch in enumerate(batchwise(
                                  screed.open(ifile, parse_description=False),
                                  batch_size)):

        norm.total += batch_size
        total = norm.total
        discarded = norm.discarded

        if index > 0 and index % 100000 == 0:
            print('... kept {kept} of {total} or {perc:2}%'
                  .format(kept=total - discarded,
                          total=total,
                          perc=int(100. - discarded / float(total) * 100.)),
                  file=sys.stderr)

            print('... in file ' + input_filename, file=sys.stderr)

            if report_fp:
                print(total + " " + total - discarded + " " +
                      1. - (discarded / float(total)), file=fp)
                report_fp.flush()

        # If in paired mode, check that the reads are properly interleaved
        if paired:
            if not check_is_pair(batch[0], batch[1]):
                raise IOError('Error: Improperly interleaved pairs \
                    {b0} {b1}'.format(b0=batch[0].name, b1=batch[1].name))

        yield batch


class Normalizer(object):
    def __init__(self, desired_coverage, htable, report_fp=None):
        self.htable = htable
        self.desired_coverage = desired_coverage
        self.report_fp = report_fp

        self.total = 0
        self.discarded = 0
        self.corrupt_files = []

    def __call__(self, input_filename, force_paired=False):
        seq = ""

        desired_coverage = self.desired_coverage
        ksize = self.htable.ksize()

        batch_size = 1
        if force_paired:
            batch_size = 2

        for batch in WithDiagnostics(input_filename, batch_size,
                                     self.report_fp, force_paired, self):
            passed_filter = False
            passed_length = True

            for record in batch:
                if len(record.sequence) < ksize:
                    passed_length = False
                    continue

                seq = record.sequence.replace('N', 'A')
                med, _, _ = self.htable.get_median_count(seq)

                if med < desired_coverage:
                    passed_filter = True

            if passed_length and passed_filter:
                for record in batch:
                    seq = record.sequence.replace('N', 'A')
                    self.htable.consume(seq)
                    yield record
            else:
                self.discarded += batch_size


def handle_error(error, output_name, input_name, fail_save, htable):
    print('** ERROR: ' + str(error), file=sys.stderr)
    print('** Failed on {name}: '.format(name=input_name), file=sys.stderr)
    if fail_save:
        tablename = os.path.basename(input_name) + '.ct.failed'
        print('** ...dumping k-mer counting table to {tn}'
              .format(tn=tablename), file=sys.stderr)
        htable.save(tablename)
    try:
        os.remove(output_name)
    except:  # pylint: disable=bare-except
        print('** ERROR: problem removing corrupt filtered file',
              file=sys.stderr)


@contextmanager
def CatchIOErrors(ifile, ofile, save_on_fail, ht, force, norm):
    """
    Context manager to do boilerplate excepting of IOErrors; also does
    upkeep on some statistics and diagnostic output.
    """

    # global corrupt_files, total, discarded, total_acc, discarded_acc
    caught_error = False
    try:
        yield
    except IOError as err:
        caught_error = True
        handle_error(err, ofile, ifile, save_on_fail, ht)
        if not force:
            print >> sys.stderr, '** Exiting!'

            sys.exit(1)
        else:
            print('*** Skipping error file, moving on...', file=sys.stderr)
            norm.corrupt_files.append(ifile)

    if norm.total == 0 and norm.discarded == 0:
        print('SKIPPED empty file ' + ifile, file=sys.stderr)
    elif not caught_error:
        total = norm.total
        discarded = norm.discarded
        print('DONE with {inp}; kept {kept} of {total} or {perc:2}%'
              .format(inp=ifile, kept=total - discarded,
                      total=total,
                      perc=int(100. - discarded / float(total) * 100.)),
              file=sys.stderr)
        print('output in ' + ofile.name, file=sys.stderr)


def normalize_by_median_and_check(input_filename, htable, single_output_file,
                                  fail_save, paired, force, norm,
                                  report_fp=None):
    total = 0
    discarded = 0

    total_acc = None
    discarded_acc = None

    if single_output_file:
        if single_output_file is sys.stdout:
            output_name = '/dev/stdout'
        else:
            output_name = single_output_file.name
        outfp = single_output_file

    else:
        output_name = os.path.basename(input_filename) + '.keep'
        outfp = open(output_name, 'w')

    with CatchIOErrors(input_filename, outfp, fail_save, htable, force, norm):

        for record in norm(input_filename, paired):
            write_record(record, outfp)

        total = norm.total
        discarded = norm.discarded

        if report_fp:
            print(str(total) + " " + str(total - discarded) + " " +
                  str(1. - (discarded / float(total))), file=report_fp)
            report_fp.flush()

    return norm.total, norm.discarded, norm.corrupt_files


def get_parser():
    epilog = ("""
    Discard sequences based on whether or not their median k-mer abundance lies
    above a specified cutoff. Kept sequences will be placed in <fileN>.keep.

    Paired end reads will be considered together if :option:`-p` is set. If
    either read will be kept, then both will be kept. This should result in
    keeping (or discarding) each sequencing fragment. This helps with retention
    of repeats, especially. With :option: `-u`/:option:`--unpaired-reads`, 
    unpaired reads from the specified file will be read after the paired data
    is read. 

    With :option:`-s`/:option:`--savetable`, the k-mer counting table
    will be saved to the specified file after all sequences have been
    processed. With :option:`-d`, the k-mer counting table will be
    saved every d files for multifile runs; if :option:`-s` is set,
    the specified name will be used, and if not, the name `backup.ct`
    will be used.  :option:`-l`/:option:`--loadtable` will load the
    specified k-mer counting table before processing the specified
    files.  Note that these tables are are in the same format as those
    produced by :program:`load-into-counting.py` and consumed by
    :program:`abundance-dist.py`.

    :option:`-f`/:option:`--fault-tolerant` will force the program to continue
    upon encountering a formatting error in a sequence file; the k-mer counting
    table up to that point will be dumped, and processing will continue on the
    next file.

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
    parser.add_argument('-p', '--paired', action='store_true')
    parser.add_argument('-u', '--unpaired-reads',
                        metavar="unpaired_reads_filename", help='with paired data only,\
                        include an unpaired file')
    parser.add_argument('-s', '--savetable', metavar="filename", default='',
                        help='save the k-mer counting table to disk after all'
                        'reads are loaded.')
    parser.add_argument('-R', '--report',
                        metavar='filename', type=argparse.FileType('w'))
    parser.add_argument('-f', '--fault-tolerant', dest='force',
                        help='continue on next file if read errors are \
                         encountered', action='store_true')
    parser.add_argument('--save-on-failure', dest='fail_save',
                        action='store_false', default=True,
                        help='Save k-mer counting table when an error occurs')
    parser.add_argument('-d', '--dump-frequency', dest='dump_frequency',
                        type=int, help='dump k-mer counting table every d '
                        'files', default=-1)
    parser.add_argument('-o', '--out', metavar="filename",
                        dest='single_output_file',
                        type=argparse.FileType('w'),
                        default=None, help='only output a single file with '
                        'the specified filename; use a single dash "-" to '
                        'specify that output should go to STDOUT (the '
                        'terminal)')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')
    parser.add_argument('--report-total-kmers', '-t', action='store_true',
                        help="Prints the total number of k-mers"
                        " post-normalization to stderr")
    parser.add_argument('--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    add_loadhash_args(parser)
    return parser


def CheckpointCountingTable(input_filenames, freq, ht, savename):
    """
    Generator/context manager to progressively save counting tables
    """

    for index, ifile in enumerate(input_filenames):
        yield ifile
        if freq > 0 and index > 0 and index % freq == 0:
            print('Backup: Saving k-mer counting file through ' +
                  ifile, file=sys.stderr)
            if savename:
                hashname = savename
                print('...saving to ' + hashname, file=sys.stderr)
            else:
                hashname = 'backup.ct'
                print('Nothing given for savetable, saving to ' + hashname,
                      file=sys.stderr)
            ht.save(hashname)


def main():  # pylint: disable=too-many-branches,too-many-statements
    info('normalize-by-median.py', ['diginorm'])
    args = get_parser().parse_args()

    report_on_config(args)

    report_fp = args.report

    # check for similar filenames
    filenames = []
    for pathfilename in args.input_filenames:
        filename = pathfilename.split('/')[-1]
        if (filename in filenames):
            print("WARNING: At least two input files are named \
%s . (The script normalize-by-median.py can not handle this, only one .keep \
file for one of the input files will be generated.)" % filename,
                  file=sys.stderr)
        else:
            filenames.append(filename)

    # check for others
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savetable:
        check_space_for_hashtable(
            args.n_tables * args.min_tablesize, args.force)

    # list to save error files along with throwing exceptions
    corrupt_files = []

    if args.loadtable:
        print('loading k-mer counting table from ' + args.loadtable,
              file=sys.stderr)
        htable = khmer.load_counting_hash(args.loadtable)
    else:
        print('making k-mer counting table', file=sys.stderr)
        htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                         args.n_tables)

    input_filename = None

    norm = Normalizer(args.cutoff, htable, report_fp)

    for f in CheckpointCountingTable(args.input_filenames, args.dump_frequency,
                                     htable, args.savetable):
        total_acc, discarded_acc, corrupt = \
            normalize_by_median_and_check(
                f, htable, args.single_output_file,
                args.fail_save, args.paired, args.force, norm, report_fp)
        corrupt_files += corrupt

    # Stuff to handle paired and unpaired data
    if args.paired and args.unpaired_reads:
        args.paired = False
        output_name = args.unpaired_reads
        if not args.single_output_file:
            output_name = os.path.basename(args.unpaired_reads) + '.keep'
        outfp = open(output_name, 'w')
        total_acc, discarded_acc, corrupt = \
            normalize_by_median_and_check(
                args.unpaired_reads, htable, args.single_output_file,
                args.fail_save, args.paired, args.force, norm, report_fp)
        corrupt_files += corrupt

    if args.report_total_kmers:
        print('Total number of unique k-mers: {0}'
              .format(htable.n_unique_kmers()),
              file=sys.stderr)

    if args.savetable:
        print('...saving to ' + args.savetable, file=sys.stderr)
        htable.save(args.savetable)

    fp_rate = \
        khmer.calc_expected_collisions(htable, args.force, max_false_pos=.8)
    # for max_false_pos see Zhang et al., http://arxiv.org/abs/1309.2975

    print('fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate),
          file=sys.stderr)

    if args.force and len(corrupt_files) > 0:
        print("** WARNING: Finished with errors!", file=sys.stderr)
        print("** IOErrors occurred in the following files:", file=sys.stderr)
        print("\t", " ".join(corrupt_files), file=sys.stderr)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
