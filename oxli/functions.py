#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#


from collections import namedtuple
import threading
import math
import khmer.utils


def optimal_size(K, M=None, f=None):
    """
    Utility function for estimating optimal counting table args where:
      - N: number of unique kmers [required]
      - M: the allotted amount of memory [optional, conflicts with f]
      - f: the desired false positive rate [optional, conflicts with M]
    """
    if all((K is not None, M is not None, f is None)):
        return estimate_optimal_with_K_and_M(K, M)
    elif all((K is not None, M is None, f is not None)):
        return estimate_optimal_with_K_and_f(K, f)
    else:
        raise TypeError("K and either M or f must be defined.")


def estimate_optimal_with_K_and_M(K, M):
    """
    Utility function for estimating optimal counting table args where K is the
    number of unique kmer and M is the allotted amount of memory
    """

    N = math.log(2) * (M / float(K))
    intN = int(N)
    if intN == 0:
        intN = 1
    X = int(M / intN)
    M = X * intN
    f2 = (1 - math.exp(-K / float(X))) ** intN
    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(intN, X, M, f2)


def estimate_optimal_with_K_and_f(K, f):
    """
    Utility function for estimating optimal memory where K is the number of
    unique kmers and f is the desired false positive rate
    """
    N = math.log(f, 0.5)
    intN = int(N)
    if intN == 0:
        intN = 1

    H1 = int(-K / (math.log(1 - f ** (1 / float(intN)))))
    M1 = H1 * intN
    f1 = (1 - math.exp(-K / float(H1))) ** intN

    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(intN, H1, M1, f1)


def optimal_args_output_gen(unique_kmers, fp_rate):
    """
    Assembles output string for optimal arg sandbox scripts
    """
    to_print = []

    to_print.append('')  # blank line
    to_print.append('number of unique k-mers: \t{0}'.format(unique_kmers))
    to_print.append('false positive rate: \t{:>.3f}'.format(fp_rate))
    to_print.append('')  # blank line
    to_print.append('If you have expected false positive rate to achieve:')
    to_print.append('expected_fp\tnumber_hashtable(Z)\tsize_hashtable(H)\t'
                    'expected_memory_usage')

    for fp_rate in range(1, 10):
        Z, H, M, f = optimal_size(unique_kmers, f=fp_rate/10.0)
        to_print.append('{:11.3f}\t{:19}\t{:17e}\t{:21e}'.format(f, Z, H, M))

    mem_list = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000, 5000]

    to_print.append('')  # blank line
    to_print.append('If you have expected memory to use:')
    to_print.append('expected_memory_usage\tnumber_hashtable(Z)\t'
                    'size_hashtable(H)\texpected_fp')

    for mem in mem_list:
        Z, H, M, f = optimal_size(unique_kmers, M=mem*1000000000)
        to_print.append('{:21e}\t{:19}\t{:17e}\t{:11.3f}'.format(M, Z, H, f))
    return "\n".join(to_print)


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
