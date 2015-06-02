#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2013-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
% python scripts/collect-variants.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.

TODO: add to sandbox README
"""
from __future__ import print_function

import sys
import screed
import os
import khmer
from khmer.khmer_args import build_counting_args

DEFAULT_NORMALIZE_LIMIT = 20


def main():
    parser = build_counting_args()
    parser.add_argument("-t", "--trusted-cutoff", dest="trusted_cutoff",
                        type=int, default=3)
    parser.add_argument("--bits-theta", help="Tuning parameter controlling"
                        "trade off of speed vs alignment sensitivity",
                        default=1.0, type=float, dest="bits_theta")
    parser.add_argument('--normalize-to', '-Z', type=int, dest='normalize_to',
                        help='base cutoff on abundance',
                        default=DEFAULT_NORMALIZE_LIMIT)
    parser.add_argument('-s', '--savehash', dest='savehash', default='')
    parser.add_argument('-l', '--loadhash', dest='loadhash',
                        default='')
    parser.add_argument('--details-out', dest="details_out")
    parser.add_argument('input_filenames', nargs='+')

    args = parser.parse_args()

    if not args.quiet:
        print('\nPARAMETERS:', file=sys.stderr)
        print(' - kmer size =    %d \t\t(-k)' % args.ksize, file=sys.stderr)
        print(' - n hashes =     %d \t\t(-N)' % args.n_tables, file=sys.stderr)
        print(' - min hashsize = %-5.2g \t(-x)' % \
            args.min_tablesize, file=sys.stderr)
        print('', file=sys.stderr)
        print('Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize)' % \
            (args.n_tables * args.min_tablesize), file=sys.stderr)
        print('-' * 8, file=sys.stderr)

    K = args.ksize
    HT_SIZE = args.min_tablesize
    N_HT = args.n_tables

    filenames = args.input_filenames

    if args.loadhash:
        print('loading hashtable from', args.loadhash)
        ht = khmer.load_counting_hash(args.loadhash)
    else:
        print('making hashtable')
        ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    aligner = khmer.ReadAligner(ht, args.trusted_cutoff, args.bits_theta)

    if args.details_out is not None:
        details_out = open(args.details_out, "w")
    else:
        details_out = None

    total = 0
    discarded = 0
    for input_filename in filenames:
        output_name = os.path.basename(input_filename) + '.keepvar'
        outfp = open(output_name, 'w')

        for n, record in enumerate(screed.open(input_filename)):
            if n > 0 and n % 10000 == 0:
                print('... kept', total - discarded, 'of', total, ', or', \
                    int(100. - discarded / float(total) * 100.), '%')
                print('... in file', input_filename)

            total += 1

            if len(record.sequence) < K:
                continue

            seq = record.sequence.upper().replace('N', 'A')

            ##

            # build the alignment...
            score, graph_alignment, read_alignment, truncated = \
                aligner.align(record.sequence)

            # next, decide whether or to keep it.
            keep = False
            if truncated:
                keep = True             # keep all truncated alignments - why?
            else:

                # build a better sequence -- this is the corrected one.
                graph_seq = graph_alignment.replace("-", "")
                # OR?
                #graph_seq = ""
                #for i in range(len(graph_alignment)):
                #    if graph_alignment[i] == "-":
                #        graph_seq += read_alignment[i]
                #    else:
                #        graph_seq += graph_alignment[i]

                # get the minimum count for this new sequence
                mincount = ht.get_min_count(graph_seq)
                if mincount < args.normalize_to:
                    keep = True

            if details_out is not None:
                details_out.write(
                    "+{7}\t{0:0.2f}\t{3}\t{4}\nread:      "
                    "{6}\ngraph_aln: {1}\nread_aln:  {2}\nstored_seq:{5}\n"
                    "".format(
                        score, graph_alignment, read_alignment, truncated,
                        keep, seq, record.sequence, record.name))

            if keep:
                ht.consume(seq)
                outfp.write('>%s\n%s\n' % (record.name, record.sequence))
            else:
                discarded += 1

        if total:
            print('DONE with', input_filename, \
                '; kept', total - discarded, 'of', total, 'or', \
                int(100. - discarded / float(total) * 100.), '%')
        print('output in', output_name)

    if args.savehash:
        print('Saving hashfile through', input_filename)
        print('...saving to', args.savehash)
        ht.save(args.savehash)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht, args.force, max_false_pos=.2)
    print('fp rate estimated to be %1.3f' % fp_rate)


if __name__ == '__main__':
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
