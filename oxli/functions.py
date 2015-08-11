"""
A collection of functions for use throughout khmer/oxli
"""

#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#


from __future__ import print_function
from collections import namedtuple
import threading
import math
import khmer.utils
import sys


def build_graph(ifilenames, graph, num_threads=1, tags=False):
    """
    Algorithm to construct a counting graph from a set of input files
    takes in list of input files, existing graph
    optionally, number of threads and if there should be tags
    """

    if tags:
        eat = graph.consume_fasta_and_tag_with_reads_parser
    else:
        eat = graph.consume_fasta_with_reads_parser

    for _, ifile in enumerate(ifilenames):
        rparser = khmer.ReadParser(ifile)
        threads = []

        for _ in range(num_threads):
            cur_thread = threading.Thread(target=eat, args=(rparser,))
            threads.append(cur_thread)
            cur_thread.start()

        for thread in threads:
            thread.join()
