#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Build a graph from the given sequences, save in <ptname>.

% python scripts/load-graph.py <ptname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys
import threading

import khmer
import argparse
from khmer.khmer_args import build_hashbits_args
from khmer.khmer_args import (report_on_config, info, add_threading_args)
from khmer.kfile import check_file_status, check_space
from khmer.kfile import check_space_for_hashtable
from oxli import build_graph, common


def get_parser():
    parser = argparse.ArgumentParser(description="Load the sequences into the "
                                     "compressible graph format plus "
                                     "optional tagset.")

    build_graph.add_args(parser)
    common.add_hashbits_args(parser)
    return parser


def main():
    info('load-graph.py', ['graph', 'SeqAn'])
    args = get_parser().parse_args()
    report_on_config(args)

    # handle argparse file things in a way that functions can handle
    ifiles = []
    for element in args.input_filenames:
        ifiles.append(element)

    build_graph.do_build_graph(args.output_filename, ifiles, args.force,
                               args.no_build_tagset, args.report_total_kmers,
                               args.write_fp_rate, args.quiet,
                               args.ksize, args.n_tables,
                               args.min_tablesize, args.threads)

if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
