#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#


import threading
import math
import khmer.utils


def estimate_optimal_with_N_and_M(N, M):
    """
    Utility function for estimating optimal counting table args where N is the
    number of unique kmer and M is the allotted amount of memory
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
    Utility function for estimating optimal memory where N is the number of
    unique kmers and f is the desired false positive rate
    """
    Z = math.log(f, 0.5)
    intZ = int(Z)
    if intZ == 0:
        intZ = 1

    H1 = int(-N/(math.log(1-f**(1/float(intZ)))))
    M1 = H1 * intZ
    f1 = (1-math.exp(-N/float(H1)))**intZ

    return intZ, H1, M1, f1


def optimal_args_output_gen(unique_kmers, fp_rate):
    """
    Assembles output string for optimal arg sandbox scripts
    """
    to_print = 'number of unique k-mers:'
    to_print += '\n' + 'false positive rate:    {:>.3f}'.format(unique_kmers,
                                                                fp_rate)
    to_print += """\n\n\nIf you have expected false positive rate to achieve:
        \nexpected_fp
        \tnumber_hashtable(Z)
        \tsize_hashtable(H)
        \texpected_memory_usage"""

    for fp_rate in range(1, 10):
        Z, H, M, f = estimate_optimal_with_N_and_f(unique_kmers, fp_rate/10.0)
        to_print = to_print + '{:11.3f}\t{:19}\t{:17e}\t{:21e}\n'.format(
            f, Z, H, M)

    mem_list = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000, 5000]

    to_print = to_print + \
        '\nIf you have expected memory to use:\n' + \
        'expected_memory_usage\tnumber_hashtable(Z)\t' + \
        'size_hashtable(H)\texpected_fp\n'
    for mem in mem_list:
        Z, H, M, f = estimate_optimal_with_N_and_M(unique_kmers,
                                                   mem*1000000000)
        to_print = to_print + '{:21e}\t{:19}\t{:17e}\t{:11.3f}\n'.format(
            M, Z, H, f)
    return to_print


def build_graph(ifilenames, graph, num_threads=1, tags=False):
    """
    Algorithm to construct a counting graph from a set of input files
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
