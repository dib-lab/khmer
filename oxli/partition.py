#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,no-member
"""Common functions for partitioning."""

import sys
import gc
import os.path

# stdlib queue module was renamed on Python 3
try:
    import queue
except ImportError:
    import Queue as queue


def worker(que, basename, stop_big_traversals):
    while True:
        try:
            (nodegraph, index, start, stop) = que.get(False)
        except queue.Empty:
            print('exiting', file=sys.stderr)
            return

        outfile = basename + '.subset.%d.pmap' % (index,)
        if os.path.exists(outfile):
            print('SKIPPING', outfile, ' -- already exists', file=sys.stderr)
            continue

        print('starting:', basename, index, file=sys.stderr)

        # pay attention to stoptags when partitioning; take command line
        # direction on whether or not to exhaustively traverse.
        subset = nodegraph.do_subset_partition(start, stop, True,
                                               stop_big_traversals)

        print('saving:', basename, index, file=sys.stderr)
        subset.save_partitionmap(outfile)
        del subset
        gc.collect()

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
