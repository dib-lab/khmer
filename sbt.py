"""
A trial implementation of sequence bloom trees, Solomon & Kingsford, 2015.

This is a simple in-memory version where all of the graphs are in
memory at once; to move it onto disk, the graphs would need to be
dynamically loaded for each query.

To try it out, do::

    factory = GraphFactory(ksize, tablesizes)
    root = Node(factory)

    graph1 = factory.create_nodegraph()
    # ... add stuff to graph1 ...
    leaf1 = Leaf("a", graph1)
    root.add_node(leaf1)

For example, ::

    # filenames: list of fa/fq files
    # ksize: k-mer size
    # tablesizes: Bloom filter table sizes

    factory = GraphFactory(ksize, tablesizes)
    root = Node(factory)

    for filename in filenames:
        graph = factory.create_nodegraph()
        graph.consume_fasta(filename)
        leaf = Leaf(filename, graph)
        root.add_node(leaf)

then define a search function, ::

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.graph.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * len(seq)):
            return 1
        return 0
"""

import random
import khmer

class GraphFactory(object):
    "Build new nodegraphs (Bloom filters) of a specific (fixed) size."

    def __init__(self, ksize, tablesizes):
        self.ksize = ksize
        self.tablesizes = tablesizes

    def create_nodegraph(self):
        return khmer._Hashbits(self.ksize, self.tablesizes)


class Node(object):
    "Internal node of SBT; has 0, 1, or 2 children."

    def __init__(self, factory):
        self.factory = factory
        self.graph = factory.create_nodegraph()
        self.children = 0
        
        self.subnodes = []

    def add_node(self, node):
        # do we have room for another child? if so, add.
        if len(self.subnodes) < 2:
            self.subnodes.append(node)
            self.graph.update(node.graph)
            self.children += 1
        # nope - insert a new node.
        else:
            if self.subnodes[0].children == self.subnodes[1].children:
                subn = random.choice(self.subnodes)
            elif self.subnodes[0].children < self.subnodes[1].children:
                subn = self.subnodes[0]
            else:
                subn = self.subnodes[1]

            ## push child down one level in the tree ##
                
            # remove from immediate:
            self.subnodes.remove(subn)

            # create new child node & fill:
            n = Node(self.factory)
            n.add_node(node)
            n.add_node(subn)

            # add new child node to ourselves
            self.subnodes.append(n)

            # don't forget to update from the new Bloom filter/nodegraph.
            self.graph.update(node.graph)
            # note: subn.graph is already included.

            # lots more children.
            self.children += 2

    def find(self, search_fn, *args):
        if not search_fn(self, *args):
            return []
        else:
            x = []
            for n in self.subnodes:
                x.extend(n.find(search_fn, *args))
            return x

class Leaf(object):
    def __init__(self, metadata, nodegraph):
        self.metadata = metadata
        self.graph = nodegraph
        self.children = 0

    def find(self, search_fn, *args):
        if search_fn(self, *args):
            return [self]
        return []


def test_simple():
    factory = GraphFactory(5, [101, 103, 117])
    root = Node(factory)

    leaf1 = Leaf("a", factory.create_nodegraph())
    leaf1.graph.count('AAAAA')
    leaf1.graph.count('AAAAT')
    leaf1.graph.count('AAAAC')
    
    leaf2 = Leaf("b", factory.create_nodegraph())
    leaf2.graph.count('AAAAA')
    leaf2.graph.count('AAAAT')
    leaf2.graph.count('AAAAG')

    leaf3 = Leaf("c", factory.create_nodegraph())
    leaf3.graph.count('AAAAA')
    leaf3.graph.count('AAAAT')
    leaf3.graph.count('CAAAA')

    leaf4 = Leaf("d", factory.create_nodegraph())
    leaf4.graph.count('AAAAA')
    leaf4.graph.count('CAAAA')
    leaf4.graph.count('GAAAA')

    leaf5 = Leaf("e", factory.create_nodegraph())
    leaf5.graph.count('AAAAA')
    leaf5.graph.count('AAAAT')
    leaf5.graph.count('GAAAA')
    
    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def search_kmer(obj, seq):
        return obj.graph.get(seq)

    leaves = [leaf1, leaf2, leaf3, leaf4, leaf5 ]
    kmers = [ "AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA" ]

    def search_kmer_in_list(kmer):
        x = []
        for l in leaves:
            if l.graph.get(kmer):
                x.append(l)

        return set(x)

    for kmer in kmers:
        assert set(root.find(search_kmer, kmer)) == search_kmer_in_list(kmer)

    print '-----'
    print [ x.metadata for x in root.find(search_kmer, "AAAAA") ]
    print [ x.metadata for x in root.find(search_kmer, "AAAAT") ]
    print [ x.metadata for x in root.find(search_kmer, "AAAAG") ]
    print [ x.metadata for x in root.find(search_kmer, "CAAAA") ]
    print [ x.metadata for x in root.find(search_kmer, "GAAAA") ]

def test_longer_search():
    ksize = 5
    factory = GraphFactory(ksize, [101, 103, 117])
    root = Node(factory)

    leaf1 = Leaf("a", factory.create_nodegraph())
    leaf1.graph.count('AAAAA')
    leaf1.graph.count('AAAAT')
    leaf1.graph.count('AAAAC')
    
    leaf2 = Leaf("b", factory.create_nodegraph())
    leaf2.graph.count('AAAAA')
    leaf2.graph.count('AAAAT')
    leaf2.graph.count('AAAAG')

    leaf3 = Leaf("c", factory.create_nodegraph())
    leaf3.graph.count('AAAAA')
    leaf3.graph.count('AAAAT')
    leaf3.graph.count('CAAAA')

    leaf4 = Leaf("d", factory.create_nodegraph())
    leaf4.graph.count('AAAAA')
    leaf4.graph.count('CAAAA')
    leaf4.graph.count('GAAAA')

    leaf5 = Leaf("e", factory.create_nodegraph())
    leaf5.graph.count('AAAAA')
    leaf5.graph.count('AAAAT')
    leaf5.graph.count('GAAAA')
    
    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.graph.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * (len(seq) - ksize + 1)):
            return 1
        return 0

    try1 = [ x.metadata for x in root.find(search_transcript, "AAAAT", 1.0) ]
    assert set(try1) == set([ 'a', 'b', 'c', 'e' ]), try1 # no 'd'
    
    try2 = [ x.metadata for x in root.find(search_transcript, "GAAAAAT", 0.6) ]
    assert set(try2) == set([ 'a', 'b', 'c', 'd', 'e' ])

    try3 = [ x.metadata for x in root.find(search_transcript, "GAAAA", 1.0) ]
    assert set(try3) == set([ 'd', 'e' ]), try3
