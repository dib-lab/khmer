#!/usr/bin/env python

import argparse
import os
from math import floor
import pickle
import sys

import khmer


class BandedHashBuffer(object):
    """
    Store k-mers in separate bands, using disk to keep memory usage low.

    The `self.hashlists` variable is the primary member. It is a dictionary
    that stores a list of hash values for each band.

    When a hash value is added to the buffer, it is stored in the list
    corresponding to the appropriate band. Once the buffer accumulates a user-
    specified number of hash values, each list is serialized and written to
    temporary files on disk, and the buffer is flushed and reset.

    Once all the input is consumed, the temporary files are re-opened and used
    to populate a k-mer count table for each band, one band at a time.
    """
    def __init__(self, numbands, outfmt='tmp.band{band}.buffer{buffer}.pickle',
                 maxsize=int(1e8)):
        self.numbands = numbands
        self.outfmt = outfmt
        self.maxsize = maxsize
        self._num_flushes = 0
        self._hash_count = 0
        self.hashlists = dict()
        self.reset()

    def __len__(self):
        # sum([len(hl) for hl in self.hashlists])
        return self._hash_count

    def reset(self):
        del(self.hashlists)
        self._hash_count = 0
        self.hashlists = dict()
        for i in range(self.numbands):
            self.hashlists[i] = list()

    def flush(self):
        if len(self) == 0:
            return

        # Keep track of how many times the buffer has been flushed to help with
        # managing temporary files.
        self._num_flushes += 1
        print('DEBUG flush', self._num_flushes, file=sys.stderr)

        for i in range(self.numbands):
            outfilename = self.outfmt.format(band=i+1,
                                             buffer=self._num_flushes)
            with open(outfilename, 'wb') as outfile:
                pickle.dump(self.hashlists[i], outfile)
        self.reset()

    def add(self, hashval):
        band = floor(hashval / (2**64) * self.numbands)
        self.hashlists[band].append(hashval)
        self._hash_count += 1
        if self._hash_count >= self.maxsize:
            self.flush()

    def get_counts(self, memory, ksize):
        if len(self) > 0:
            self.flush()

        for band in range(self.numbands):
            counts = khmer.Counttable(ksize, memory / 4, 4)
            for i in range(self._num_flushes):
                bufferfilename = self.outfmt.format(band=band+1, buffer=i+1)
                with open(bufferfilename, 'rb') as bufferfile:
                    hashlist = pickle.load(bufferfile)
                    for khash in hashlist:
                        counts.add(khash)
                os.remove(bufferfilename)
            yield band, counts


parser = argparse.ArgumentParser()
parser.add_argument('-k', '--ksize', type=int, metavar='K', default=31,
                    help='k-mer size')
parser.add_argument('-n', '--num-bands', type=int, metavar='N',
                    help='number of bands')
parser.add_argument('-b', '--buffersize', type=float, metavar='B', default=1e8,
                    help='number of k-mers to keep in memory before writing '
                    'the buffer to disk and flushing')
parser.add_argument('-m', '--memory', type=float, metavar='M', default=1e4,
                    help='memory (in bytes) to allocate to each counttable in '
                    'the output')
parser.add_argument('-o', '--outfmt', metavar='FMT', default='band{}.ct',
                    help='a string template for output files; default is '
                    '"band{}.ct"; brackets will be replaced with band numbers')
parser.add_argument('infiles', nargs='+')
args = parser.parse_args()


kg = khmer.Counttable(args.ksize, 1, 1)
kbuffer = BandedHashBuffer(args.num_bands, maxsize=int(args.buffersize))
for infile in args.infiles:
    reads = khmer.ReadParser(infile)
    for read in reads:
        for kmer in kg.get_kmer_hashes(read.sequence):
            kbuffer.add(kmer)

for band, counttable in kbuffer.get_counts(args.memory, args.ksize):
    fpr = khmer.calc_expected_collisions(counttable, max_false_pos=100.0)
    ctfilename = args.outfmt.format(band + 1)
    print('Band', band, 'FPR', fpr, ctfilename, file=sys.stderr)
    counttable.save(ctfilename)
