#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Estimate optimal arguments using HLL couter.

% python sandbox/optimal_args_HLL.py [ -k <k size> ] [ -e <ERROR_RATE> ] <data1> <data2> ...
Use '-h' for parameter help.
"""


import argparse
import os
import sys
import textwrap
import math
import khmer
from khmer.khmer_args import DEFAULT_K, info, ComboFormatter
from khmer import __version__


def get_parser():
    descr = "Estimate optimal arguments using HLL couter., with precision <= ERROR_RATE."
    epilog = ("""
    A HyperLogLog counter is used to do cardinality estimation. Since this counter
    is based on a tradeoff between precision and memory consumption,
    :option:`-e`/:option:`--error-rate` can be used to control how much
    memory will be used. In practice the memory footprint is small even
    at low error rates (< 0.01).
    :option:`-k`/:option:`--ksize` should be set to the desired k-mer size.
    Output is sent to STDOUT, but a report file can be generated with
    :option:`-R`/:option:`--report`.
    Example::
        optimal_args_HLL.py -k 17 tests/test-data/test-abund-read{,-2,-3}.fa
    Example::
""" "    optimal_args_HLL.py.py -R optimal_args -k 30 tests/test-data/test-abund-read-paired.fa")  # noqa
    parser = argparse.ArgumentParser(
        description=descr, epilog=textwrap.dedent(epilog),
        formatter_class=ComboFormatter)

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)

    parser.add_argument('--version', action='version',
                        version='khmer {v}'.format(v=__version__))
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')

    parser.add_argument('--ksize', '-k', type=int, default=env_ksize,
                        help='k-mer size to use')

    parser.add_argument('--error-rate', '-e', type=float, default=0.01,
                        help='Acceptable error rate')

    parser.add_argument('-R', '--report',
                        metavar='filename', type=argparse.FileType('w'))

    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')


    return parser

def estimate_optimal_with_N_and_M(N,M):
    Z = math.log(2)*(M/float(N))
    intZ = int(Z)
    if intZ == 0:
        intZ = 1
    H = int(M/intZ)
    M = H*intZ
    f2 = (1-math.exp(-N/float(H)))**intZ
    return intZ,H,M,f2

def estimate_optimal_with_N_and_f(N,f):
    Z = math.log(f,0.5)
    intZ = int(Z)
    if intZ == 0:
        intZ = 1
    H1 = int(-N/(math.log(1-f**(1/float(intZ)))))
    M1 = H1 * intZ
    f1 = (1-math.exp(-N/float(H1)))**intZ
    return intZ, H1,M1,f1

def to_print_func(unique_kmers,fp_rate):
    to_print = 'number of unique k-mers:    {}\n\
    false positive rate:    {:>.3f}'.format(unique_kmers,fp_rate) 
    to_print = to_print + \
    '\n\n\nIf you have expected false positive rate to achieve:\n' + \
'expected_fp\tnumber_hashtable(Z)\tsize_hashtable(H)\texpected_memory_usage\n'

    for fp_rate in range(1,10):
        Z,H,M,f = estimate_optimal_with_N_and_f(unique_kmers,fp_rate/10.0)
        to_print =to_print + '{:11.3f}\t{:19}\t{:17e}\t{:21e}\n'.format(f,Z,H,M)

    mem_list = [1,5,10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000,5000]
    
    to_print = to_print + \
               '\nIf you have expected memory to use:\n' + \
                'expected_memory_usage\tnumber_hashtable(Z)\t' + \
                'size_hashtable(H)\texpected_fp\n'
    for mem in mem_list:
        Z,H,M,f = estimate_optimal_with_N_and_M(unique_kmers,mem*1000000000)
        #print Z,H,M,f
        to_print =to_print + '{:21e}\t{:19}\t{:17e}\t{:11.3f}\n'.format(M,Z,H,f) 
    return to_print

def main():
    info('optimal_args_HLL.py', ['SeqAn'])
    args = get_parser().parse_args()

    hllcpp = khmer.HLLCounter(args.error_rate, args.ksize)

    report_fp = args.report
    fp_rate = args.error_rate
    input_filename = None
    for index, input_filename in enumerate(args.input_filenames):
        hllcpp.consume_fasta(input_filename)

    cardinality = hllcpp.estimate_cardinality()
    print >> sys.stdout, 'Estimated number of unique k-mers: {0}'.format(
        cardinality)
        
    to_print = to_print_func(cardinality,fp_rate)
    
    print >> sys.stdout, to_print
    if report_fp:
        print >> report_fp, to_print
        report_fp.flush()
    
    



if __name__ == "__main__":
    main()