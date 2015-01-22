#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2010-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Produce the k-mer abundance distribution for the given file, without
loading a prebuilt k-mer counting table.

% python scripts/abundance-dist-single.py <data> <histout>

Use '-h' for parameter help.
"""
import os
import sys
import khmer
import threading
import textwrap
from khmer.khmer_args import (build_counting_args, add_threading_args,
                              report_on_config, info)
from oxli import common, abund_dist_single
import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="Caluculate the abundance "
                                     "distribution of k-mers from a single "
                                     "sequence file",
                                     epilog=textwrap.dedent(
                                         abund_dist_single.parser_epilog),
                                     formatter_class=common.ComboFormatter)
    # use modified build_counting_args and build_hash_args from oxli.common
    common.add_counting_args(parser)
    abund_dist_single.add_args(parser)

    return parser


def main():  # pylint: disable=too-many-locals,too-many-branches
    info('abundance-dist-single.py', ['counting', 'SeqAn'])
    args = get_parser().parse_args()
    report_on_config(args)

    abund_dist_single.do_abund_dist_single(args.input_sequence_filename,
                                           args.output_histogram_filename,
                                           args.output_zero, args.bigcount,
                                           args.squash_output, args.savetable,
                                           args.report_total_kmers,
                                           args.quiet, args.ksize,
                                           args.n_tables,
                                           args.min_tablesize,
                                           args.threads, args.force)



if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
