#! /usr/bin/env python
"""
Build an outline index for the given sequence file.

% python scripts/build-outline-index.py @CTB

Use '-h' for parameter help.
"""

import khmer
import argparse

DEFAULT_K = 32

def main():
    parser = argparse.ArgumentParser(
        description="Build an outline index for sequence retrieval via k-mer")

    parser.add_argument('graphname')
    parser.add_argument('seqdb')

    parser.add_argument('-k', '--ksize', dest='ksize', default=DEFAULT_K)
    
    args = parser.parse_args()
    graphname = args.graphname
    seqdb = args.seqdb
    K = args.ksize
    
    # define empty hashtable -- we don't need the hashtable for building the
    # outline index, only for arbitrary k-mer queries.
    ht = khmer.new_hashbits(K, 1, 1)

    # load the tagset
    ht.load_tagset(graphname + '.tagset')

    # now, binary index the file
    bin_filename = khmer.convert_fasta_to_indexed_bin(seqdb)
    
    # build the outline index for the file using the tags
    ht.build_outline_index(bin_filename)

if __name__ == '__main__':
    main()
