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
The output is in JSON format by default, simplifying parsing/piping for third party programs and pipelines, ie:

{"sample": ["sample.fastq"], "filter": "eschColi_K12.pt", "contamination": 0.3}

Usage:
% python scripts/filter-contamination.py <htname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
from __future__ import division

# Built the test graph using:
# scripts/load-graph.py -x 4e9 -N 4 --no-build-tagset foo.khmer ../data/16s.fa

import sys

import khmer
import json
from collections import defaultdict
import argparse

from khmer.khmer_args import info
from khmer.file import check_file_status, check_space


def get_parser():
    parser = argparse.ArgumentParser(
        description="Determine which organisms are present in a given file ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'graph', help="basename of the input k-mer presence table")
    parser.add_argument('data', help="files to be decontaminated")
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser

# http://scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/


def sliding_window_it(sequence, winSize, step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""

    # Verify the inputs
    try:
        it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(0) == type(winSize)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception(
            "**ERROR** winSize must not be larger than sequence length.")

    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence) - winSize) / step) + 1

    # Do the work
    for i in range(0, int(numOfChunks * step), step):
        yield sequence[i:i + winSize]


def main():
    info('filter-contamination.py', ['graph', 'classification'])
    args = get_parser().parse_args()
    graph = args.graph
    data = args.data

    results = defaultdict(dict)

    # Get samples to be analyzed against the graph file
    filenames = [data]
    for _ in filenames:
        check_file_status(_)

    # Is there enough disk space available for input/output files before doing
    # anything?
    check_space(filenames)

    print 'loading ht graph %s' % graph
    htable = khmer.load_hashbits(graph)

    for _, filename in enumerate(filenames):
        print('querying sample {sample} against filter {filt}').format(
            sample=filename, filt=graph)
        ksize = htable.ksize()
        total_query_kmers = 0
        contaminant_total_matches = 0

        rparser = khmer.ReadParser(filename, 1)

        for r in rparser:
            read_kmers = len(r.sequence) - ksize + 1
            contaminant_read_matches = 0

            for position in range(read_kmers):
                for kmer in sliding_window_it(r.sequence, ksize):
                    contaminant_read_matches += htable.get(kmer)

            contaminant_total_matches += contaminant_read_matches
            total_query_kmers += read_kmers

        contam = contaminant_read_matches / total_query_kmers

        results['sample'] = filenames
        results['filter'] = graph
        results['contamination'] = contam

        print json.dumps(results)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
