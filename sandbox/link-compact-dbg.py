#! /usr/bin/env python
import khmer
import screed
import argparse
from collections import OrderedDict
import sys

# graph settings
DEFAULT_KSIZE=31
NODEGRAPH_SIZE=8e8 # small, big is 2e8

# minhash settings
MH_SIZE_DIVISOR=50
MH_MIN_SIZE=5
MH_PRIME=9999999967

class Pathfinder(object):
    "Track segment IDs & adjacency lists."
    def __init__(self, ksize):
        self.ksize = ksize

        self.segment_counter = 1
        self.segments = {}
        self.segments_r = {}
        self.adjacencies = {}

    def new_segment(self, kmer):
        if kmer in self.segments_r:
            return self.segments_r[kmer]

        this_id = self.segment_counter
        self.segment_counter += 1

        self.segments[this_id] = self.ksize
        self.segments_r[kmer] = this_id

        return this_id

    def new_linear_segment(self, size):
        this_id = self.segment_counter
        self.segment_counter += 1
        self.segments[this_id] = size
        return this_id

    def add_adjacency(self, node_id, adj):
        node_id, adj = min(node_id, adj), max(node_id, adj)

        x = self.adjacencies.get(node_id, set())
        x.add(adj)
        self.adjacencies[node_id] = x


def traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy, degree_nodes,
                                   lh):
    size, conns, visited = graph.traverse_linear_path(nk, degree_nodes,
                                                      stop_bf)
    if not size:
        return

    linear_path_labels = set()
    for node in visited:
        linear_path_labels.update(lh.get_tag_labels(node))

    # TODO: do something graph-y with these, like split paths across HDNs
    # or simply omit HDNs altogether with traversal.

    # give it a segment ID
    path_id = pathy.new_linear_segment(size)

    # for all adjacencies, add.
    for conn in conns:
        conn_id = pathy.segments_r.get(conn)
        pathy.add_adjacency(path_id, conn_id)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfiles', nargs='+')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    parser.add_argument('-x', '--tablesize', default=NODEGRAPH_SIZE,
                            type=float)
    parser.add_argument('--force', action='store_true')
    #parser.add_argument('--gml', action='store_true')
    args = parser.parse_args()

    assert args.ksize % 2, "ksize must be odd"
    assert args.output, "you probably want an output file"

    print('building graphs and loading files')

    # Create graph, and two stop bloom filters - one for loading, one for
    # traversing. Create them all here so that we can error out quickly
    # if memory is a problem.

    graph = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    stop_bf = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    stop_bf2 = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    n = 0

    # load in all of the input sequences, one file at a time.
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            n += 1
            if n % 10000 == 0:
                print('...', seqfile, n)
            graph.consume(record.sequence)

    # complain if too small set of graphs was used.
    fp_rate = khmer.calc_expected_collisions(graph,
                                             args.force, max_false_pos=.05)

    # initialize the object that will track information for us.
    pathy = Pathfinder(args.ksize)

    print('finding high degree nodes')
    degree_nodes = khmer.HashSet(args.ksize)
    n = 0
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            n += 1
            if n % 10000 == 0:
                print('...2', seqfile, n)
            # walk across sequences, find all high degree nodes,
            # name them and cherish them. Don't do this on identical sequences.
            if min(stop_bf2.get_kmer_counts(record.sequence)) == 0:
                stop_bf2.consume(record.sequence)
                degree_nodes += graph.find_high_degree_nodes(record.sequence)
    del stop_bf2

    if not len(degree_nodes):
        print('no high degree nodes; exiting.')
        sys.exit(0)

    ####

    lh = khmer._GraphLabels(graph)
    n = 0
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            n += 1
            if n % 10000 == 0:
                print('...2', seqfile, n)
            lh.label_across_high_degree_nodes(record.sequence, degree_nodes, n)

    print('num labels:', lh.n_labels())

    # get all of the degree > 2 nodes and give them IDs.
    for node in degree_nodes:
        pathy.new_segment(node)

    print('traversing linear segments from', len(degree_nodes), 'nodes')

    # now traverse from each high degree nodes into all neighboring nodes,
    # seeking adjacencies.  if neighbor is high degree node, add it to
    # adjacencies; if neighbor is not, then traverse the linear path.  also
    # track minhashes while we're at it.
    for n, k in enumerate(degree_nodes):
        if n % 10000 == 0:
            print('...', n, 'of', len(degree_nodes))

        # retrieve the segment ID of the primary node.
        k_id = pathy.segments_r[k]

        # find all the neighbors of this high-degree node.
        nbh = graph.neighbors(k)
        for nk in nbh:
            # neighbor is high degree? fine, mark its adjacencies.
            if nk in degree_nodes:
                nk_id = pathy.segments_r[nk]
                pathy.add_adjacency(k_id, nk_id)
            else:
                # linear! walk it.
                traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy,
                                               degree_nodes, lh)

    print(len(pathy.segments), 'segments, containing',
              sum(pathy.segments.values()), 'nodes')

    # save to GML
    if args.output:
        import graph_writer

        print('saving to', args.output)
        fp = open(args.output, 'w')
        w = graph_writer.GmlWriter(fp, [], [])

        for k, v in pathy.segments.items():
            w.add_vertex(k, v, [])

        for k, v in pathy.adjacencies.items():
            for edge in v:
                w.add_edge(k, edge, [])


if __name__ == '__main__':
    main()
