#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2010-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
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
Extract partitioned sequences into files grouped by partition size.

% python scripts/extract-partitions.py <base> <file1.part> [ <file2.part> ... ]

Grouped sequences will be <base>.groupN.fa files.

Use '-h' for parameter help.

@CTB note that if threshold is != 1, those sequences will not be output
by output_unassigned...
"""

import sys
import screed
import textwrap
from contextlib import contextmanager
import khmer

from khmer.kfile import (check_input_files, check_space,
                         add_output_compression_type,
                         get_file_writer)
from khmer.khmer_args import sanitize_help, KhmerArgumentParser
from khmer.utils import write_record

DEFAULT_MAX_SIZE = int(1e6)
DEFAULT_THRESHOLD = 5


def read_partition_file(filename):
    """Utility function to get partitioned reads from file."""
    for record_index, record in enumerate(screed.open(filename)):
        _, partition_id = record.name.rsplit('\t', 1)
        yield record_index, record, int(partition_id)


def get_parser():
    """Create parser for extract-partitions.py."""
    epilog = """
    Example (results will be in ``example.group0000.fa``)::

        load-graph.py -k 20 example tests/test-data/random-20-a.fa
        partition-graph.py example
        merge-partitions.py -k 20 example
        annotate-partitions.py -k 20 example tests/test-data/random-20-a.fa
        extract-partitions.py example random-20-a.fa.part

    (:program:`extract-partitions.py` will produce a partition size
    distribution in <base>.dist. The columns are: (1) number of reads,
    (2) count of partitions with n reads, (3) cumulative sum of partitions,
    (4) cumulative sum of reads.)
    """
    parser = KhmerArgumentParser(
        description="Separate sequences that are annotated with partitions "
        "into grouped files.", epilog=textwrap.dedent(epilog),
        citations=['graph'])
    parser.add_argument('prefix', metavar='output_filename_prefix')
    parser.add_argument('part_filenames', metavar='input_partition_filename',
                        nargs='+')
    parser.add_argument('-X', '--max-size', dest='max_size',
                        default=DEFAULT_MAX_SIZE, type=int,
                        help='Max group size (n sequences)')
    parser.add_argument('-m', '--min-partition-size', dest='min_part_size',
                        default=DEFAULT_THRESHOLD, type=int,
                        help='Minimum partition size worth keeping')
    parser.add_argument('-n', '--no-output-groups', dest='output_groups',
                        default=True, action='store_false',
                        help='Do not actually output groups files.')
    parser.add_argument('-U', '--output-unassigned', default=False,
                        action='store_true',
                        help='Output unassigned sequences, too')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    add_output_compression_type(parser)
    return parser


@contextmanager
def PartitionedReader(filename_list, quiet=False, single=False):
    yield PartitionedReadIterator(filename_list, quiet, single)


def PartitionedReadIterator(filename_list, quiet=False, single=False):
    """
    Generator to do boilerplate output of statistics.

    Uses a list of input files and verbosity
    Returns reads and partition IDs
    """
    for filename in filename_list:
        for index, read, pid in read_partition_file(filename):
            if not quiet:
                if index % 100000 == 0:
                    print('...x2', index, file=sys.stderr)
            yield read, pid
            if single:
                break  # only yield a single read from each file


class PartitionExtractor(object):
    """
    Does extraction, processing and accounting of partitioned reads.

    Contains methods for self_maintenance  and output production
    """

    def __init__(self, file_list, min_size, max_size):
        # We'll make our own generator! With files! And...
        self.file_list = file_list
        self.n_unassigned = 0
        self.count = {}

        self.divvy = None
        self.group_n = 0
        self.group_d = {}
        self.min_size = min_size
        self.max_size = max_size

    def process_unassigned(self, outfp=None):
        """
        Process unassigned reads.

        Can optionally output said reads if outfp is given
        Also develops counts of partition IDs--necessary for further processing
        """
        with PartitionedReader(self.file_list) as reader:
            for read, pid in reader:
                self.count[pid] = self.count.get(pid, 0) + 1

                if pid == 0:
                    self.n_unassigned += 1
                    if outfp:
                        write_record(read, outfp)

    def output_histogram(self, dist_filename):
        """Output histogram of partition counts to the given filename."""
        # develop histogram of partition sizes
        dist = {}
        for _, size in list(self.count.items()):
            dist[size] = dist.get(size, 0) + 1

        # output histogram
        distfp = open(dist_filename, 'w')

        total = 0
        wtotal = 0
        for counter, index in sorted(dist.items()):
            total += index
            wtotal += counter * index
            distfp.write('%d %d %d %d\n' % (counter, index, total, wtotal))
        distfp.close()

    def develop_groups(self):
        """Processing method that divides up the partitions into groups."""
        if 0 in self.count:            # eliminate unpartitioned sequences
            del self.count[0]

        # sort groups by size
        self.divvy = sorted(list(self.count.items()), key=lambda y: y[1])
        self.divvy = [y for y in self.divvy if y[1] > self.min_size]

        # divvy up into different groups, based on having max_size sequences
        # in each group.
        total = 0
        group = set()
        for partition_id, n_reads in self.divvy:
            group.add(partition_id)
            total += n_reads

            if total > self.max_size:
                for partition_id in group:
                    self.group_d[partition_id] = self.group_n

                self.group_n += 1
                group = set()
                total = 0

        if group:
            for partition_id in group:
                self.group_d[partition_id] = self.group_n
            self.group_n += 1

    class ReadGroupGenerator(object):
        """
        Generator that yields partitioned reads and their group.

        Takes PartitionExtractor and PartitionedReadIterator objects
        """

        def __init__(self, extractor):
            self.extractor = extractor

            self.total_seqs = 0
            self.part_seqs = 0
            self.toosmall_parts = 0

        def __call__(self, reader):
            for read, partition_id in reader:
                self.total_seqs += 1
                if partition_id == 0:
                    continue

                try:
                    group_n = self.extractor.group_d[partition_id]
                except KeyError:
                    assert self.extractor.count[partition_id] <=\
                        self.extractor.min_size
                    self.toosmall_parts += 1
                    continue

                yield read, group_n
                self.part_seqs += 1


def main():
    args = sanitize_help(get_parser()).parse_args()

    distfilename = args.prefix + '.dist'

    for infile in args.part_filenames:
        check_input_files(infile, args.force)

    check_space(args.part_filenames, args.force)

    print('---', file=sys.stderr)
    print('reading partitioned files:', repr(
        args.part_filenames), file=sys.stderr)
    if args.output_groups:
        print('outputting to files named "%s.groupN.fa"' %
              args.prefix, file=sys.stderr)
        print('min reads to keep a partition:',
              args.min_part_size, file=sys.stderr)
        print('max size of a group file:', args.max_size, file=sys.stderr)
    else:
        print('NOT outputting groups! Beware!', file=sys.stderr)

    if args.output_unassigned:
        print('outputting unassigned reads to "%s.unassigned.fa"' %
              args.prefix, file=sys.stderr)
    print('partition size distribution will go to %s'
          % distfilename, file=sys.stderr)
    print('---', file=sys.stderr)

    #

    suffix = None
    is_fastq = None

    with PartitionedReader(args.part_filenames, True, True) as reader:
        for read, _ in reader:
            if is_fastq is None:
                is_fastq = hasattr(read, 'quality')
            else:
                assert hasattr(read, 'quality') == is_fastq,\
                    "Input files must have consistent format."

    if is_fastq:
        suffix = "fq"
    else:
        suffix = "fa"

    # remember folks, generators exhaust themseleves
    extractor = PartitionExtractor(args.part_filenames,
                                   args.min_part_size,
                                   args.max_size)

    if args.output_unassigned:
        ofile = open('%s.unassigned.%s' % (args.prefix, suffix), 'wb')
        unassigned_fp = get_file_writer(ofile, args.gzip, args.bzip)
        extractor.process_unassigned(unassigned_fp)
        unassigned_fp.close()
    else:
        extractor.process_unassigned()

    extractor.output_histogram(distfilename)

    if not args.output_groups:
        sys.exit(0)

    extractor.develop_groups()

    print('%d groups' % extractor.group_n, file=sys.stderr)
    if extractor.group_n == 0:
        print('nothing to output; exiting!', file=sys.stderr)
        return

    # open a bunch of output files for the different groups
    group_fps = {}
    for index in range(extractor.group_n):
        fname = '%s.group%04d.%s' % (args.prefix, index, suffix)
        group_fp = get_file_writer(open(fname, 'wb'), args.gzip,
                                   args.bzip)
        group_fps[index] = group_fp

    # write 'em all out!
    # refresh the generator
    read_generator = PartitionExtractor.ReadGroupGenerator(extractor)

    with PartitionedReader(args.part_filenames) as reader:
        for read, group_n in read_generator(reader):
            outfp = group_fps[group_n]
            write_record(read, outfp)

    print('---', file=sys.stderr)
    print('Of %d total seqs,' % read_generator.total_seqs, file=sys.stderr)
    print('extracted %d partitioned seqs into group files,' %
          read_generator.part_seqs, file=sys.stderr)
    print('discarded %d sequences from small partitions (see -m),' %
          read_generator.toosmall_parts, file=sys.stderr)
    print('and found %d unpartitioned sequences (see -U).' %
          extractor.n_unassigned, file=sys.stderr)
    print('', file=sys.stderr)
    print('Created %d group files named %s.groupXXXX.%s' %
          (len(group_fps),
           args.prefix,
           suffix), file=sys.stderr)


if __name__ == '__main__':
    main()
