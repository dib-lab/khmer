#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring


import threading
import math
import khmer.utils


def estimate_optimal_with_N_and_M(N, M):
    """
    Utility function for estimating optimal memory
    """
    Z = math.log(2)*(M/float(N))
    intZ = int(Z)
    if intZ == 0:
        intZ = 1
    H = int(M/intZ)
    M = H*intZ
    f2 = (1-math.exp(-N/float(H)))**intZ
    return intZ, H, M, f2


def estimate_optimal_with_N_and_f(N, f):
    """
    Utility function for estimating optimal memory
    """
    Z = math.log(f, 0.5)
    intZ = int(Z)
    if intZ == 0:
        intZ = 1

    H1 = int(-N/(math.log(1-f**(1/float(intZ)))))
    M1 = H1 * intZ
    f1 = (1-math.exp(-N/float(H1)))**intZ

    return intZ, H1, M1, f1


def build_graph(ifilenames, graph, num_threads=1, tags=False):

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
