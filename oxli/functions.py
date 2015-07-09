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


def estimate_optimal_with_N_and_M(N, M):
    """
    Utility function for estimating optimal counting table args where N is the
    number of unique kmer and M is the allotted amount of memory
    """
    Z = math.log(2) * (M / float(N))
    intZ = int(Z)
    if intZ == 0:
        intZ = 1
    H = int(M / intZ)
    M = H * intZ
    f2 = (1 - math.exp(-N / float(H))) ** intZ
    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(intZ, H, M, f2)


def estimate_optimal_with_N_and_f(N, f):
    """
    Utility function for estimating optimal memory where N is the number of
    unique kmers and f is the desired false positive rate
    """
    Z = math.log(f, 0.5)
    intZ = int(Z)
    if intZ == 0:
        intZ = 1

    H1 = int(-N / (math.log(1 - f ** (1 / float(intZ)))))
    M1 = H1 * intZ
    f1 = (1 - math.exp(-N / float(H1))) ** intZ

    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(intZ, H1, M1, f1)


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
        Z, H, M, f = estimate_optimal_with_N_and_f(
            unique_kmers, fp_rate / 10.0)
        to_print.append('{:11.3f}\t{:19}\t{:17e}\t{:21e}'.format(f, Z, H, M))

    mem_list = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000, 5000]

    to_print.append('')  # blank line
    to_print.append('If you have expected memory to use:')
    to_print.append('expected_memory_usage\tnumber_hashtable(Z)\t'
                    'size_hashtable(H)\texpected_fp')

    for mem in mem_list:
        Z, H, M, f = estimate_optimal_with_N_and_M(unique_kmers,
                                                   mem * 1000000000)
        to_print.append('{:21e}\t{:19}\t{:17e}\t{:11.3f}'.format(M, Z, H, f))
    return "\n".join(to_print)


def do_sanity_checking(args, desired_max_fp):
    # if optimization args are given, do optimization
    if args.unique_kmers != 0:
        if args.max_memory_usage:
            # verify that this is a sane memory usage restriction
            res = estimate_optimal_with_N_and_M(args.unique_kmers,
                                                args.max_memory_usage)
            if res.fp_rate > desired_max_fp:
                print("""
*** ERROR: The given restrictions yield an estimate false positive rate of {0},
*** which is above the recommended false positive ceiling of 0.1!
*** Aborting!""".format(res.fp_rate), file=sys.stderr)
                sys.exit(1)
        else:
            res = estimate_optimal_with_N_and_f(args.unique_kmers,
                                                desired_max_fp)
            if args.max_tablesize and args.max_tablesize < res.htable_size:
                print("*** Warning: The given tablesize is be too small!",
                      file=sys.stderr)
                print("*** Estimating false positive rate to be {0}".format(
                      res.fp_rate), file=sys.stderr)
            else:
                print("*** INFO: set memory ceiling using atuo optimaztion.",
                      file=sys.stderr)
                print("*** Ceiling is: {0} bytes\n".format(res.mem_use),
                      file=sys.stderr)
                args.max_mem = res.mem_use

    return args


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
