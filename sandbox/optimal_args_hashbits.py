#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=invalid-name,missing-docstring
"""
Estimate optimal arguments using hashbits counting.

% python sandbox/optimal_args_hashbits.py  <data1> [ <data2> <...> ]

Use '-h' for parameter help.
"""
from __future__ import print_function

import sys
import math
import threading

import khmer
from khmer.khmer_args import build_hashbits_args
from khmer.khmer_args import (report_on_config, info, add_threading_args)
from khmer.kfile import check_input_files, check_space
from khmer.kfile import check_space


def get_parser():
    parser = build_hashbits_args(descr="Load sequences into the compressible "
                                 "graph format plus optional tagset.")
    add_threading_args(parser)
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        nargs='+', help='input FAST[AQ] sequence filename')
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
    info('optimal_args_hashbits.py', ['graph', 'SeqAn'])
    args = get_parser().parse_args()
    report_on_config(args, hashtype='hashbits')


    filenames = args.input_filenames
    base = filenames[0]
    for _ in args.input_filenames:
        check_input_files(_, False)

    check_space(args.input_filenames, False)
    check_space((float(args.n_tables * args.min_tablesize) / 8.), False)

    print('Counting kmers from sequences in %s' % repr(filenames),
          file=sys.stderr)

    htable = khmer.new_hashbits(args.ksize, args.min_tablesize, args.n_tables)
    target_method = htable.consume_fasta_with_reads_parser

    for _, filename in enumerate(filenames):
        rparser = khmer.ReadParser(filename)
        threads = []
        print('consuming input', filename, file=sys.stderr)
        for num in xrange(args.threads):
            cur_thread = threading.Thread(
                target=target_method, args=(rparser,))
            threads.append(cur_thread)
            cur_thread.start()

        for thread in threads:
            thread.join()
    unique_kmers = htable.n_unique_kmers()
    print('Total number of unique k-mers: {0}'.format(unique_kmers),
          file=sys.stderr)

    info_optimal = open(base + '.optimal_args', 'w')

    fp_rate = khmer.calc_expected_collisions(htable)
    print('fp rate estimated to be %1.3f' % fp_rate, file=sys.stderr)

    if fp_rate > 0.15:          # 0.18 is ACTUAL MAX. Do not change.
        print("**", file=sys.stderr)
        print("** ERROR: the graph structure is too small for this data set."
              "Increase table size/# tables.", file=sys.stderr)
        print("**", file=sys.stderr)
        if not False:
            sys.exit(1)

    to_print = to_print_func(unique_kmers,fp_rate)
    
    print(to_print, file=info_optimal)
    
    print('optimal arguments were written to', base + '.optimal_args',
          file=sys.stderr)
    
if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
