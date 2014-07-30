#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Screen for contamination in a similar way as FACS, fastq_screen or deconseq do.

% python scripts/filter-contamination.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

# Built the test graph using:
# scripts/load-graph.py -x 4e9 -N 4 --no-build-tagset foo.khmer ../data/16s.fa

import sys

import khmer
import argparse

from khmer.khmer_args import info
from khmer.file import check_file_status, check_space

def get_parser():
    parser = argparse.ArgumentParser(
        description="Determine which organisms are present in a given file ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('graph', help="basename of the input k-mer presence table")
    parser.add_argument('data', help="files to be decontaminated")
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    info('filter-contamination.py', ['graph', 'classification'])
    args = get_parser().parse_args()
    graph = args.graph
    data = args.data

    # Get samples to be analyzed against the graph file
    filenames = [data]
    for _ in filenames:
        check_file_status(_)

    # Is there enough disk space available for input/output files before doing anything?
    check_space(filenames)

    print 'loading ht graph %s' % graph
    htable = khmer.load_hashbits(graph)

    n_reads = 0
    len_reads = []

    # For each read determine % k-mers that are present in the hashbits table
    # XXX: Which number of **reads** from input (data argument) are present in the hashbits table (graph argument)?
    #
    # Output should resemble:
    # {
    # "input": "200_sequenced_reads.fastq"
    # "filter": "eschColi_K12_reference_genome.pt"
    # "total_read_count": 201,
    # "contaminated_reads": 1,
    # "contamination_rate": 0.004975,
    #}

    for _, filename in enumerate(filenames):
        rparser = khmer.ReadParser(filename, 1)
        print 'consuming input {0}\n'.format(filename)

        for r in rparser:
            len_reads.append(len(r.sequence))
            n_reads = n_reads + 1
            htable.consume(r.sequence)

    print >> sys.stderr, 'Total number reads: {0}'.format(n_reads)
    print >> sys.stderr, 'Average read length: {0}'.format(sum(len_reads)/n_reads)
    print >> sys.stderr, 'Total number of occupied k-mers: {0}'.format(htable.n_occupied())
    print >> sys.stderr, 'Total number of unique k-mers: {0}'.format(htable.n_unique_kmers())

    # XXX: Make a decision based on a threshold and write to disk in appropriate places

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
