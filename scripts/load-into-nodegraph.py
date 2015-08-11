#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
# pylint disable=missing-docstring
"""
Build a graph from the given sequences, save in <ptname>.

% python scripts/load-into-nodegraph.py <ptname> <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""

import sys

from khmer.khmer_args import build_hashbits_args
from oxli import build_graph


def get_parser():
    parser = build_hashbits_args(descr="Load sequences into the compressible "
                                       "graph format plus optional tagset.")

    parser = build_graph.build_parser(parser)
    return parser


if __name__ == '__main__':
    build_graph.main(get_parser().parse_args())
    sys.exit(0)

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
