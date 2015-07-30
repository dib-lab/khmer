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


def optimal_size(num_kmers, mem_cap=None, fp_rate=None):
    """
    Utility function for estimating optimal counting table args where:
      - num_kmers: number of unique kmers [required]
      - mem_cap: the allotted amount of memory [optional, conflicts with f]
      - fp_rate: the desired false positive rate [optional, conflicts with M]
    """
    if all((num_kmers is not None, mem_cap is not None, fp_rate is None)):
        return estimate_optimal_with_K_and_M(num_kmers, mem_cap)
    elif all((num_kmers is not None, mem_cap is None, fp_rate is not None)):
        return estimate_optimal_with_K_and_f(num_kmers, fp_rate)
    else:
        raise TypeError("num_kmers and either mem_cap or fp_rate"
                        " must be defined.")


def estimate_optimal_with_K_and_M(num_kmers, mem_cap):
    """
    Utility function for estimating optimal counting table args where num_kmers
    is the number of unique kmer and mem_cap is the allotted amount of memory
    """

    n_tables = math.log(2) * (mem_cap / float(num_kmers))
    int_n_tables = int(n_tables)
    if int_n_tables == 0:
        int_n_tables = 1
    ht_size = int(mem_cap / int_n_tables)
    mem_cap = ht_size * int_n_tables
    fp_rate = (1 - math.exp(-num_kmers / float(ht_size))) ** int_n_tables
    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(int_n_tables, ht_size, mem_cap, fp_rate)


def estimate_optimal_with_K_and_f(num_kmers, des_fp_rate):
    """
    Utility function for estimating optimal memory where num_kmers  is the
    number of unique kmers and des_fp_rate is the desired false positive rate
    """
    n_tables = math.log(des_fp_rate, 0.5)
    int_n_tables = int(n_tables)
    if int_n_tables == 0:
        int_n_tables = 1

    ht_size = int(-num_kmers / (
        math.log(1 - des_fp_rate ** (1 / float(int_n_tables)))))
    mem_cap = ht_size * int_n_tables
    fp_rate = (1 - math.exp(-num_kmers / float(ht_size))) ** int_n_tables

    res = namedtuple("result", ["num_htables", "htable_size", "mem_use",
                                "fp_rate"])
    return res(int_n_tables, ht_size, mem_cap, fp_rate)


def optimal_args_output_gen(unique_kmers, fp_rate):
    """
    Assembles output string for optimal arg sandbox scripts
    takes in unique_kmers and desired fp_rate
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
        num_tables, table_size, mem_cap, fp_rate = \
            optimal_size(unique_kmers, fp_rate=fp_rate / 10.0)
        to_print.append('{:11.3f}\t{:19}\t{:17e}\t{:21e}'.format(fp_rate,
                                                                 num_tables,
                                                                 table_size,
                                                                 mem_cap))

    mem_list = [1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000, 5000]

    to_print.append('')  # blank line
    to_print.append('If you have expected memory to use:')
    to_print.append('expected_memory_usage\tnumber_hashtable(Z)\t'
                    'size_hashtable(H)\texpected_fp')

    for mem in mem_list:
        num_tables, table_size, mem_cap, fp_rate =\
            optimal_size(unique_kmers, mem_cap=mem * 1000000000)
        to_print.append('{:21e}\t{:19}\t{:17e}\t{:11.3f}'.format(mem_cap,
                                                                 num_tables,
                                                                 table_size,
                                                                 fp_rate))
    return "\n".join(to_print)


def check_fp_rate(args, desired_max_fp):
    """
    simple function to check if the restrictions in the args (if there are any)
    make sense--If not, complain. If no restrictions are given, add some that
    make sense.
    Takes in args and desired max FP rate
    """
    if args.unique_kmers != 0:
        if args.max_memory_usage:
            # verify that this is a sane memory usage restriction
            res = estimate_optimal_with_K_and_M(args.unique_kmers,
                                                args.max_memory_usage)
            if res.fp_rate > desired_max_fp:
                print("""
*** ERROR: The given restrictions yield an estimate false positive rate of {0},
*** which is above the recommended false positive ceiling of {1}!"""
                      .format(res.fp_rate, desired_max_fp), file=sys.stderr)
                if not args.force:
                    print("NOTE: This can be overridden using the --force"
                          " argument", file=sys.stderr)
                    print("*** Aborting...!", file=sys.stderr)
                    sys.exit(1)
        else:
            res = estimate_optimal_with_K_and_f(args.unique_kmers,
                                                desired_max_fp)
            if args.max_tablesize and args.max_tablesize < res.htable_size:
                print("*** Warning: The given tablesize is too small!",
                      file=sys.stderr)
                print("*** Estimating false positive rate to be {0}".format(
                    res.fp_rate), file=sys.stderr)
            else:
                print("*** INFO: set memory ceiling automatically.",
                      file=sys.stderr)
                print("*** Ceiling is: {0} bytes\n".format(res.mem_use),
                      file=sys.stderr)
                args.max_mem = res.mem_use

    return args


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
