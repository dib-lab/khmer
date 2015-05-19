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
import screed
from khmer.khmer_args import build_hashbits_args
from khmer.khmer_args import (report_on_config, info, add_threading_args)
from khmer.kfile import check_input_files, check_space
from khmer.kfile import check_space_for_hashtable
from oxli import khmer_api, build_graph

def get_parser():
    parser = build_hashbits_args(descr="Load-Graph called.")
    parser = build_graph.build_parser(parser)
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    build_graph.main(args)

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
