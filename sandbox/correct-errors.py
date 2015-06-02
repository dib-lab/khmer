#! /usr/bin/env python
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Streaming error correction.

% python sandbox/correct-errors.py [ <data1> [ <data2> [ ... ] ] ]

Use -h for parameter help.

TODO: paired support: paired reads should be kept together.
TODO: load/save counting table.
TODO: move output_single elsewhere
TODO: add to sandbox/README
TODO: change name to correct-reads?
"""
from __future__ import print_function
import sys
import screed
import os
import khmer
import argparse
import tempfile
import shutil

DEFAULT_NORMALIZE_LIMIT = 20
DEFAULT_CUTOFF = 2

DEFAULT_K = 32
DEFAULT_N_HT = 4
DEFAULT_MIN_HASHSIZE = 1e6


def output_single(read, new_sequence):
    name = read.name
    sequence = new_sequence

    quality = None
    if hasattr(read, 'quality'):
        quality = read.quality[:len(sequence)]

        # in cases where sequence _lengthened_, need to truncate it to
        # match the quality score length.
        sequence = sequence[:len(quality)]

    if quality:
        assert len(sequence) == len(quality), (sequence, quality)
        return "@%s\n%s\n+\n%s\n" % (name, sequence, quality)
    else:
        return ">%s\n%s\n" % (name, sequence)


def main():
    parser = argparse.ArgumentParser(description='XXX')

    env_ksize = os.environ.get('KHMER_KSIZE', DEFAULT_K)
    env_n_hashes = os.environ.get('KHMER_N_HASHES', DEFAULT_N_HT)
    env_hashsize = os.environ.get('KHMER_MIN_HASHSIZE', DEFAULT_MIN_HASHSIZE)

    parser.add_argument('--ksize', '-k', type=int, dest='ksize',
                        default=env_ksize,
                        help='k-mer size to use')
    parser.add_argument('--n_hashes', '-N', type=int, dest='n_hashes',
                        default=env_n_hashes,
                        help='number of hash tables to use')
    parser.add_argument('--hashsize', '-x', type=float, dest='min_hashsize',
                        default=env_hashsize,
                        help='lower bound on hashsize to use')

    parser.add_argument("--trusted-cov", dest="trusted_cov", type=int,
                        default=DEFAULT_CUTOFF)
    parser.add_argument("--theta", dest="bits_theta", type=float, default=1.0)

    parser.add_argument('--normalize-to', '-Z', type=int, dest='normalize_to',
                        help='base cutoff on median k-mer abundance of this',
                        default=DEFAULT_NORMALIZE_LIMIT)

    parser.add_argument('--tempdir', '-T', type=str, dest='tempdir',
                        default='./')

    parser.add_argument('input_filenames', nargs='+')
    args = parser.parse_args()

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    NORMALIZE_LIMIT = args.normalize_to

    print('making hashtable')
    ht = khmer.new_counting_hash(K, HT_SIZE, N_HT)

    aligner = khmer.ReadAligner(ht, args.trusted_cov, args.bits_theta)

    tempdir = tempfile.mkdtemp('khmer', 'tmp', args.tempdir)
    print('created temporary directory %s; use -T to change location' % tempdir)

    ###

    save_pass2 = 0
    n_aligned = 0
    n_corrected = 0
    total_reads = 0

    pass2list = []
    for filename in args.input_filenames:
        pass2filename = os.path.basename(filename) + '.pass2'
        pass2filename = os.path.join(tempdir, pass2filename)
        corrfilename = os.path.basename(filename) + '.corr'

        pass2list.append((filename, pass2filename, corrfilename))

        pass2fp = open(pass2filename, 'w')
        corrfp = open(corrfilename, 'w')

        for n, read in enumerate(screed.open(filename)):
            total_reads += 1

            if n % 10000 == 0:
                print('...', n, filename, n_aligned, n_corrected, save_pass2, \
                      total_reads)
            seq = read.sequence.replace('N', 'A')

            # build the alignment...
            score, graph_alignment, read_alignment, truncated = \
                aligner.align(read.sequence)

            # next, decide whether or to keep it.
            output_corrected = False
            if not truncated:
                n_aligned += 1

                # build a better sequence -- this is the corrected one.
                if True:
                    graph_seq = graph_alignment.replace("-", "")
                else:
                    graph_seq = ""
                    for i in range(len(graph_alignment)):
                        if graph_alignment[i] == "-":
                            graph_seq += read_alignment[i]
                        else:
                            graph_seq += graph_alignment[i]

                corrected = graph_seq
                if graph_seq != read.sequence:
                    n_corrected += 1

                # get the minimum count for this new sequence
                mincount = ht.get_min_count(graph_seq)
                if mincount < args.normalize_to:
                    output_corrected = True

            # has this portion of the graph saturated? if not,
            # consume & save => pass2.
            if output_corrected:
                corrfp.write(output_single(read, corrected))
            else:  # uncorrected...
                ht.consume(read.sequence)
                pass2fp.write(output_single(read, read.sequence))
                save_pass2 += 1

        pass2fp.close()
        corrfp.close()

        print('%s: kept aside %d of %d from first pass, in %s' % \
              (filename, save_pass2, n, filename))
        print('aligned %d of %d reads so far' % (n_aligned, total_reads))
        print('changed %d of %d reads so far' % (n_corrected, total_reads))

    for orig_filename, pass2filename, corrfilename in pass2list:
        print('second pass: looking at sequences kept aside in %s' % \
              pass2filename)
        for n, read in enumerate(screed.open(pass2filename)):
            if n % 10000 == 0:
                print('... x 2', n, pass2filename, n_aligned, n_corrected, \
                      total_reads)

            corrfp = open(corrfilename, 'a')

            # build the alignment...
            score, graph_alignment, read_alignment, truncated = \
                aligner.align(read.sequence)

            if truncated:               # no good alignment; output original
                corrected = read.sequence
            else:
                n_aligned += 1
                # build a better sequence -- this is the corrected one.
                if True:
                    graph_seq = graph_alignment.replace("-", "")
                else:
                    graph_seq = ""
                    for i in range(len(graph_alignment)):
                        if graph_alignment[i] == "-":
                            graph_seq += read_alignment[i]
                        else:
                            graph_seq += graph_alignment[i]

                corrected = graph_seq
                if corrected != read.sequence:
                    n_corrected += 1

            corrfp.write(output_single(read, corrected))

        print('removing %s' % pass2filename)
        os.unlink(pass2filename)

    print('removing temp directory & contents (%s)' % tempdir)
    shutil.rmtree(tempdir)

    print('Aligned %d of %d total' % (n_aligned, total_reads))
    print('Changed %d of %d total' % (n_corrected, total_reads))

if __name__ == '__main__':
    main()
