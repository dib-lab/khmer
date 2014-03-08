#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Count the overlap k-mers, which are the k-mers apperaring in two sequence
datasets.

usage: count-overlap_cpp.py [-h] [-q] [--ksize KSIZE] [--n_hashes N_HASHES]
                        [--hashsize HASHSIZE]
                           1st_dataset(htfile) 2nd_dataset(fastafile) result

Use '-h' for parameter help.

"""
import khmer
from khmer.file_api import check_file_status, check_space
from khmer.khmer_args import build_hashbits_args, report_on_config
#
DEFAULT_K = 32
DEFAULT_N_HT = 4
DEFAULT_HASHSIZE = 1e6


def main():
    parser = build_hashbits_args()
    parser.add_argument('htfile')
    parser.add_argument('fafile')
    parser.add_argument('report_filename')
    args = parser.parse_args()

    report_on_config(args, hashtype='hashbits')

    K = args.ksize
    HT_SIZE = args.hashsize
    N_HT = args.n_hashes
    htfile = args.htfile
    fafile = args.fafile
    output_filename = args.report_filename
    curve_filename = output_filename + '.curve'

    infiles = [htfile, fafile]
    for infile in infiles:
        check_file_status(infile)

    check_space(infiles)

    print 'loading hashbits from', htfile
    ht1 = khmer.load_hashbits(htfile)
    K = ht1.ksize()

    output = open(output_filename, 'w')
    f_curve_obj = open(curve_filename, 'w')

    ht2 = khmer.new_hashbits(K, HT_SIZE, N_HT)

    (n_unique, n_overlap, list) = ht2.count_overlap(fafile, ht1)

    printout1 = """\
dataset1(ht file): %s
dataset2: %s

# of unique k-mers in dataset2: %d
# of overlap unique k-mers: %d

""" % (htfile, fafile, n_unique, n_overlap)
    output.write(printout1)

    figure_list1 = []
    figure_list2 = []

    for i in range(100):
        to_print = str(list[100 + i]) + ' ' + str(list[i]) + '\n'
        f_curve_obj.write(to_print)


if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
