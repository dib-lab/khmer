#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Eliminate reads with median k-mer abundance higher than
DESIRED_COVERAGE.  Output sequences will be placed in 'infile.keep'.

% python scripts/normalize-by-median.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.
"""

import sys
import screed
import os
import khmer
import textwrap
from itertools import izip
from khmer.khmer_args import (build_counting_args, add_loadhash_args,
                              report_on_config, info)
import argparse
from khmer.file import (check_space, check_space_for_hashtable,
                        check_valid_file_exists)
DEFAULT_DESIRED_COVERAGE = 10

MAX_FALSE_POSITIVE_RATE=0.8             # see Zhang et al.,
                                        # http://arxiv.org/abs/1309.2975

# Iterate a collection in arbitrary batches
# from: http://stackoverflow.com/questions/4628290/pairs-from-single-list


def batchwise(coll, size):
    iter_coll = iter(coll)
    return izip(*[iter_coll] * size)

# Returns true if the pair of records are properly pairs


def validpair(read0, read1):
    return read0.name[-1] == "1" and \
        read1.name[-1] == "2" and \
        read0.name[0:-1] == read1.name[0:-1]


# pylint: disable=too-many-locals,too-many-branches
def normalize_by_median(input_filename, outfp, htable, args, report_fp=None):

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
        if index > 0 and index % 100000 == 0:
            print '... kept {kept} of {total} or {perc:2}%'.format(
                kept=total - discarded, total=total,
                perc=int(100. - discarded / float(total) * 100.))
            print '... in file', input_filename

            if report_fp:
                print >> report_fp, total, total - discarded, \
                    1. - (discarded / float(total))
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
                htable.consume(seq)
                passed_filter = True

        # Emit records if any passed
        if passed_length and passed_filter:
            for record in batch:
                if hasattr(record, 'accuracy'):
                    outfp.write(
                        '@{name}\n{seq}\n'
                        '+\n{acc}\n'.format(name=record.name,
                                            seq=record.sequence,
                                            acc=record.accuracy))
                else:
                    outfp.write(
                        '>{name}\n{seq}\n'.format(name=record.name,
                                                  seq=record.sequence))
        else:
            discarded += batch_size

    return total, discarded


def handle_error(error, output_name, input_name, fail_save, htable):
    print >> sys.stderr, '** ERROR:', error
    print >> sys.stderr, '** Failed on {name}: '.format(name=input_name)
    if fail_save:
        tablename = os.path.basename(input_name) + '.ct.failed'
        print >> sys.stderr, \
            '** ...dumping k-mer counting table to {tn}'.format(tn=tablename)
        htable.save(tablename)
    try:
        os.remove(output_name)
    except:  # pylint: disable=bare-except
        print >> sys.stderr, '** ERROR: problem removing corrupt filtered file'


def get_parser():
    epilog = ("""
    Discard sequences based on whether or not their median k-mer abundance lies
    above a specified cutoff. Kept sequences will be placed in <fileN>.keep.

    Paired end reads will be considered together if :option:`-p` is set. If
    either read will be kept, then both will be kept. This should result in
    keeping (or discarding) each sequencing fragment. This helps with retention
    of repeats, especially.

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

    Example::

        normalize-by-median.py -k 17 tests/test-data/test-abund-read-2.fa

    Example::

""" "        normalize-by-median.py -p -k 17 tests/test-data/test-abund-read-paired.fa"  # noqa
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
    parser.add_argument('-s', '--savetable', metavar="filename", default='')
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
                        dest='single_output_filename',
                        default='', help='only output a single'
                        ' file with the specified filename')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')
    add_loadhash_args(parser)
    return parser


def main():  # pylint: disable=too-many-branches,too-many-statements
    info('normalize-by-median.py', ['diginorm'])
    args = get_parser().parse_args()

    report_on_config(args)

    report_fp = args.report

    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames)
    if args.savetable:
        check_space_for_hashtable(args.n_tables * args.min_tablesize)

    # list to save error files along with throwing exceptions
    if args.force:
        corrupt_files = []

    if args.loadtable:
        print 'loading k-mer counting table from', args.loadtable
        htable = khmer.load_counting_hash(args.loadtable)
    else:
        print 'making k-mer counting table'
        htable = khmer.new_counting_hash(args.ksize, args.min_tablesize,
                                         args.n_tables)

    total = 0
    discarded = 0

    for index, input_filename in enumerate(args.input_filenames):
        if args.single_output_filename != '':
            output_name = args.single_output_filename
            outfp = open(args.single_output_filename, 'a')
        else:
            output_name = os.path.basename(input_filename) + '.keep'
            outfp = open(output_name, 'w')

        total_acc = 0
        discarded_acc = 0

        try:
            total_acc, discarded_acc = normalize_by_median(input_filename,
                                                           outfp, htable, args,
                                                           report_fp)
        except IOError as err:
            handle_error(err, output_name, input_filename, args.fail_save,
                         htable)
            if not args.force:
                print >> sys.stderr, '** Exiting!'
                sys.exit(1)
            else:
                print >> sys.stderr, '*** Skipping error file, moving on...'
                corrupt_files.append(input_filename)
        else:
            if total_acc == 0 and discarded_acc == 0:
                print 'SKIPPED empty file', input_filename
            else:
                total += total_acc
                discarded += discarded_acc
                print 'DONE with {inp}; kept {kept} of {total} or {perc:2}%'\
                    .format(inp=input_filename,
                            kept=total - discarded, total=total,
                            perc=int(100. - discarded / float(total) * 100.))
                print 'output in', output_name

        if (args.dump_frequency > 0 and
                index > 0 and index % args.dump_frequency == 0):
            print 'Backup: Saving k-mer counting file through', input_filename
            if args.savetable:
                hashname = args.savetable
                print '...saving to', hashname
            else:
                hashname = 'backup.ct'
                print 'Nothing given for savetable, saving to', hashname
            htable.save(hashname)

    if args.savetable:
        print 'Saving k-mer counting table through', input_filename
        print '...saving to', args.savetable
        htable.save(args.savetable)

    fp_rate = khmer.calc_expected_collisions(htable)
    print 'fp rate estimated to be {fpr:1.3f}'.format(fpr=fp_rate)

    if args.force and len(corrupt_files) > 0:
        print >> sys.stderr, "** WARNING: Finished with errors!"
        print >> sys.stderr, "** IOErrors occurred in the following files:"
        print >> sys.stderr, "\t", " ".join(corrupt_files)

    if fp_rate > MAX_FALSE_POSITIVE_RATE:
        print >> sys.stderr, "**"
        print >> sys.stderr, ("** ERROR: the k-mer counting table is too small"
                              " for this data set.  Increase tablesize/# "
                              "tables.")
        print >> sys.stderr, "**"
        print >> sys.stderr, "** Do not use these results!!"
        sys.exit(1)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
