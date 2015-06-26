#! /usr/bin/env python
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
from contextlib import contextmanager

from khmer.khmer_args import (build_counting_args, add_loadhash_args,
                              report_on_config, info)
import argparse
from khmer.kfile import (check_space, check_space_for_hashtable,
                         check_valid_file_exists)
from khmer.utils import write_record, check_is_pair, broken_paired_reader


DEFAULT_DESIRED_COVERAGE = 20


def WithDiagnostics(ifilename, norm, reader, fp):
    """
    Generator/context manager to do boilerplate output of statistics using a
    Normalizer object.
    """

    index = 0

    # per read diagnostic output
    for index, record in enumerate(norm(reader)):

        if norm.total % 100000 == 0:
            print('... kept {kept} of {total} or {perc:2}% so far'
                  .format(kept=norm.total - norm.discarded,
                          total=norm.total,
                          perc=int(100. - norm.discarded /
                                   float(norm.total) * 100.)),
                  file=sys.stderr)

            print('... in file ' + ifilename, file=sys.stderr)

        yield record

    # per file diagnostic output
    if norm.total == 0:
        print('SKIPPED empty file ' + ifilename, file=sys.stderr)
    else:
        print('DONE with {inp}; kept {kept} of {total} or {perc:2}%'
              .format(inp=ifilename, kept=norm.total - norm.discarded,
                      total=norm.total, perc=int(100. - norm.discarded /
                                                 float(norm.total) * 100.)),
              file=sys.stderr)

    if fp:
        print("{total} {kept} {discarded}"
              .format(total=norm.total, kept=norm.total - norm.discarded,
                      discarded=1. - (norm.discarded / float(norm.total))),
              file=fp)
        fp.flush()


class Normalizer(object):
    """
    Digital normalization algorithm.
    """

    def __init__(self, desired_coverage, htable):
        self.htable = htable
        self.desired_coverage = desired_coverage

        self.total = 0
        self.discarded = 0

    def __call__(self, reader):
        """
        Actually does digital normalization - the core algorithm.

        * get one (unpaired) or two (paired) reads;
        * sanitize the sequences (convert Ns to As);
        * get the median k-mer count of one/both reads;
        * if any read's median k-mer count is below desired coverage, keep all;
        * consume and yield kept reads.
        """

        desired_coverage = self.desired_coverage

        for index, is_paired, read0, read1 in reader:
            passed_filter = False

            self.total += 1

            if is_paired:
                self.total += 1

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
            else:
                self.discarded += len(batch)


@contextmanager
def CatchIOErrors(ifile, out, single_out, force, corrupt_files):
    """
    Context manager to do boilerplate handling of IOErrors.
    """
    try:
        yield
    except (IOError, ValueError) as error:
        print('** ERROR: ' + str(error), file=sys.stderr)
        print('** Failed on {name}: '.format(name=ifile), file=sys.stderr)
        if not single_out:
            os.remove(out.name)
        if not force:
            print('** Exiting!', file=sys.stderr)

            sys.exit(1)
        else:
            print('*** Skipping error file, moving on...', file=sys.stderr)
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
    processed. With :option:`-d`, the k-mer counting table will be
    saved every d files for multifile runs; if :option:`-s` is set,
    the specified name will be used, and if not, the name `backup.ct`
    will be used.  :option:`-l`/:option:`--loadtable` will load the
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
                        metavar='filename', type=argparse.FileType('w'))
    parser.add_argument('-f', '--force', dest='force',
                        help='continue on next file if read errors are \
                         encountered', action='store_true')
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
    info('normalize-by-median.py', ['diginorm'])
    args = get_parser().parse_args()

    report_on_config(args)

    report_fp = args.report
    force_single = args.force_single

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
            print('ERROR: Duplicate filename--Cannot handle this!',
                  file=sys.stderr)
            print('** Exiting!', file=sys.stderr)
            sys.exit(1)

        basenames.append(basename)

    # check that files exist and there is sufficient output disk space.
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savetable:
        check_space_for_hashtable(
            args.n_tables * args.min_tablesize, args.force)

    # load or create counting table.
    if args.loadtable:
        print('loading k-mer counting table from ' + args.loadtable,
              file=sys.stderr)
        htable = khmer.load_counting_hash(args.loadtable)
    else:
        print('making k-mer counting table', file=sys.stderr)
        htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                         args.n_tables)

    input_filename = None

    # create an object to handle diginorm of all files
    norm = Normalizer(args.cutoff, htable)

    # make a list of all filenames and if they're paired or not;
    # if we don't know if they're paired, default to allowing but not
    # forcing pairing.
    files = []
    for e in filenames:
        files.append([e, args.paired])
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

    #
    # main loop: iterate over all files given, do diginorm.
    #

    for filename, require_paired in files:
        if not args.single_output_file:
            output_name = os.path.basename(filename) + '.keep'
            outfp = open(output_name, 'w')

        # failsafe context manager in case an input file breaks
        with CatchIOErrors(filename, outfp, args.single_output_file,
                           args.force, corrupt_files):

            screed_iter = screed.open(filename, parse_description=False)
            reader = broken_paired_reader(screed_iter, min_length=args.ksize,
                                          force_single=force_single,
                                          require_paired=require_paired)

            # actually do diginorm
            for record in WithDiagnostics(filename, norm, reader, report_fp):
                if record is not None:
                    write_record(record, outfp)

            print('output in ' + output_name, file=sys.stderr)
            if output_name is not '/dev/stdout':
                outfp.close()

    # finished - print out some diagnostics.

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
