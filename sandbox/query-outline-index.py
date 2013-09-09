#! /usr/bin/env python
"""
Build an outline index for the given sequence file.

% python scripts/build-outline-index.py @CTB

Use '-h' for parameter help.
"""

import sys
import khmer
import argparse
import os.path

DEFAULT_K = 32

def load_query_line(line, k):
    line = line.strip().upper()
    try:
        kmer = long(line)
    except ValueError:
        kmer = khmer.forward_hash(line, k)

    return kmer

def main():
    parser = argparse.ArgumentParser(
        description="Query the outline index to retrieve reads by k-mer.")

    parser.add_argument('graphname')
    parser.add_argument('seqdb')
    parser.add_argument('queryfile')
    parser.add_argument('--no-search-graph', '-n', default=False,
                        action='store_true', dest='no_search_graph',
                        help='Only retrieve sequences by exact tags')

    args = parser.parse_args()
    graphname = args.graphname
    seqdb = args.seqdb
    queryfile = args.queryfile

    # load the graph
    ht_filename = graphname + '.ht'
    sys.stderr.write("Loading graph/hashtable from %s\n" % ht_filename)
    ht = khmer.load_hashbits(ht_filename)

    sys.stderr.write("The K parameter is %d\n" % ht.ksize())

    # load the tagset
    tagset_filename = graphname + '.tagset'
    sys.stderr.write("Loading tagset from %s\n" % tagset_filename)
    ht.load_tagset(tagset_filename)

    # retrieve the outline index filename
    bin_filename = seqdb + '.bin'
    assert os.path.exists(bin_filename), \
        "bin filename %s doesn't exist" % (bin_filename,)

    index_filename = seqdb + '.bin.index'
    assert os.path.exists(index_filename), \
        "index filename %s doesn't exist" % (index_filename,)

    # let's load the query file.  Here, we will allow either DNA strings
    # OR hash numbers.
    
    query_set = set()

    fp = open(queryfile)
    for line in fp:
        kmer = load_query_line(line, ht.ksize())
        print 'kmer is', kmer
        query_set.add(kmer)

    sys.stderr.write("Loaded %d query k-mers\n" % len(query_set))

    tagset = None
    if args.no_search_graph:
        tagset = query_set
    else: # search the graph!
        tagset = set()
        for kmer in query_set:
            taglist = ht.find_all_tags_to_taglist(kmer)
            tagset.update(taglist)

    sys.stderr.write("Found %d tags corresponding to %d query k-mers\n" % \
                         (len(tagset), len(query_set)))
    
    read_id_list = khmer.outline_retrieve_read_ids_by_taglist(index_filename,
                                                              list(tagset))

    sys.stderr.write("Got %d read IDs\n" % len(read_id_list))

    for read_id in read_id_list:
        print 'READ ID IS', read_id
        read = khmer.outline_retrieve_read_by_id(bin_filename, read_id)
        print read
    
if __name__ == '__main__':
    main()
