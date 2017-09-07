#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name
"""
Estimate optimal choice of hash table parameters

First scenario: we know the number of k-mers N and the size of memory 
available to use M. We want to know the optimal number of hash table Z 
to get the lowest false positive rate.

For this scenario, use only "-M" and "number of kmers".
% sandbox/estimate_optimal_hash.py <number_of_kmers> -M <size_of_memory>

Second scenario: we know the number of k-mers N and the desired maximum 
false positive rate f. We want to know the minimum memory usage required 
to achieve f. 

For this scenario, use only "-f" and "number of kmers".
% sandbox/estimate_optimal_hash.py <number_of_kmers> -f <desired_fpr>

Use '-h' for parameter help.

"""
import argparse
import khmer, oxli
from khmer.khmer_args import info, optimal_size, sanitize_help
import textwrap
import sys

def get_parser():
    epilog = """

    First scenario: we know the number of k-mers N and the size of memory 
    available to use M. We want to know the optimal number of hash table Z 
    to get the lowest false positive rate.

    For this scenario, use only "-M" and "number of kmers".
    % sandbox/estimate_optimal_hash.py <number_of_kmers> -M <size_of_memory>

    Second scenario: we know the number of k-mers N and the desired maximum 
    false positive rate f. We want to know the minimum memory usage required 
    to achieve f. 

    For this scenario, use only "-f" and "number of kmers".
    % sandbox/estimate_optimal_hash.py <number_of_kmers> -f <desired_fpr>

    """
    parser = argparse.ArgumentParser(
        description='Estimate optimal choice of hash table parameters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(epilog))

    parser.add_argument('-N', help='number of estimated distinct k-mers',
                        type = int)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-M', help='size of memory available to use', 
                        type = int)
    group.add_argument('-f', help='desired maximum false posotive rate',
                        type = float)
    parser.add_argument('--version', action='version', version='%(prog)s '
                        + khmer.__version__)
    return parser


def main():
    info('estimate_optimal_hash.py', ['counting'])
    args = sanitize_help(get_parser()).parse_args()
    N = args.N
    if args.M:
        M = args.M
        result = optimal_size(N, M=M)
        print("number of estimated distinct k-mers:  ", N, file=sys.stderr)
        print("size of memory available to use:      ", M, file=sys.stderr)
        print("optimal number of hash tables:        ", result.num_htables,
              file=sys.stderr)
        print("optimal size of hash tables:          ", result.htable_size,
              file=sys.stderr)
        print("estimated false positive rate:        ", result.fp_rate,
              file=sys.stderr)
        print("estimated usage of memory:            ", result.mem_use,
              file=sys.stderr)

    elif args.f:
        f = args.f
        result = optimal_size(N, f=f)
        print("number of estimated distinct k-mers:  ", N, file=sys.stderr)
        print("desired maximum false positive rate:  ", f, file=sys.stderr)
        print("optimal number of hash tables:        ", result.num_htables,
              file=sys.stderr)
        print("optimal size of hash tables:          ", result.htable_size,
              file=sys.stderr)
        print("estimated false positive rate:        ", result.fp_rate,
              file=sys.stderr)
        print("estimated usage of memory:            ", result.mem_use,
              file=sys.stderr)
        
    else:
        get_parser().error('No action requested, add -M (size of memory available to use) or -f (desired maximum false posotive rate)')

if __name__ == '__main__':
    main()
