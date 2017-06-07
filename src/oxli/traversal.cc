/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"
#include "oxli/traversal.hh"
#include "oxli/alphabets.hh"
#include "oxli/kmer_hash.hh"

using namespace std;

namespace oxli
{

/******************************************
 * NodeGatherer
 ******************************************/

template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashgraph * ht,
                                      KmerFilterList filters) :
    KmerFactory(ht->ksize()), graph(ht), filters(filters)
{
    bitmask = 0;
    for (unsigned int i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    rc_left_shift = _ksize * 2 - 2;
}


template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashgraph * ht) :
    NodeGatherer(ht, KmerFilterList())
{
}


template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashgraph * ht,
                                      KmerFilter filter) :
    NodeGatherer(ht, KmerFilterList())
{
    filters.push_back(filter);
}


template <bool direction>
WordLength NodeGatherer<direction>::ksize() const
{
    return graph->ksize();
}


template<>
Kmer NodeGatherer<TRAVERSAL_LEFT>::get_neighbor(const Kmer& node, const char ch)
const
{
    // optimized bit-foo to check for LEFT neighbors in both forward and
    // reverse-complemented directions
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return build_kmer(kmer_f, kmer_r);
}


template<>
Kmer NodeGatherer<TRAVERSAL_RIGHT>::get_neighbor(const Kmer& node,
                                       const char ch)
const
{
    // optimized bit-foo to check for LEFT neighbors in both forward and
    // reverse-complemented directions
    HashIntoType kmer_f, kmer_r;
    kmer_f = (((node.kmer_f) << 2) & bitmask) | (twobit_repr(ch));
    kmer_r = ((node.kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift);
    return build_kmer(kmer_f, kmer_r);
}


template<bool direction>
unsigned int NodeGatherer<direction>::neighbors(const Kmer& node,
        KmerQueue & node_q)
const
{
    unsigned int found = 0;

    for (auto base : alphabets::DNA_SIMPLE) {
        // Get the putative neighboring Kmer
        Kmer neighbor = get_neighbor(node, base);
        // Now check if it's in the graph and passes the filters
        if (graph->get_count(neighbor) && !(apply_kmer_filters(neighbor, filters))) {
            node_q.push(neighbor);
            ++found;
        }
        ++base;
    }

    return found;
}


template<bool direction>
unsigned int NodeGatherer<direction>::degree(const Kmer& node)
const
{
    unsigned int degree = 0;

    for (auto base : alphabets::DNA_SIMPLE) {
        if (graph->get_count(get_neighbor(node, base))) {
            ++degree;
        }
        ++base;
    }

    return degree;
}

/******************************************
 * NodeCursor
 ******************************************/

template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashgraph * ht,
                                  Kmer start_kmer,
                                  KmerFilterList filters) :
    NodeGatherer<direction>(ht, filters)
{
    cursor = start_kmer;
}


template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashgraph * ht,
                                  Kmer start_kmer) :
    NodeCursor<direction>(ht, start_kmer, KmerFilterList())
{
}


template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashgraph * ht,
                                  Kmer start_kmer,
                                  KmerFilter filter) :
    NodeCursor<direction>(ht, start_kmer)
{
    push_filter(filter);
}


template<bool direction>
unsigned int NodeCursor<direction>::cursor_degree()
const
{
    return this->degree(this->cursor);
}



/******************************************
 * Traverser
 ******************************************/

Traverser::Traverser(const Hashgraph * ht,
                     KmerFilterList filters) :
    KmerFactory(ht->ksize()),
    graph(ht),
    left_gatherer(ht, filters),
    right_gatherer(ht, filters)
{
}

Traverser::Traverser(const Hashgraph * ht,
                     KmerFilter filter) :
    KmerFactory(ht->ksize()),
    graph(ht),
    left_gatherer(ht, filter),
    right_gatherer(ht, filter)
{
}


void Traverser::push_filter(KmerFilter filter)
{
    left_gatherer.push_filter(filter);
    right_gatherer.push_filter(filter);
}


KmerFilter Traverser::pop_filter()
{
    left_gatherer.pop_filter();
    return right_gatherer.pop_filter();
}


unsigned int Traverser::traverse(const Kmer& node,
                                 KmerQueue& node_q) const
{
    return left_gatherer.neighbors(node, node_q) +
           right_gatherer.neighbors(node, node_q);
}


unsigned int Traverser::traverse_left(const Kmer& node,
                                      KmerQueue& node_q) const
{
    return left_gatherer.neighbors(node, node_q);
}


unsigned int Traverser::traverse_right(const Kmer& node,
                                       KmerQueue& node_q) const
{
    return right_gatherer.neighbors(node, node_q);
}


unsigned int Traverser::degree(const Kmer& node) const
{
    return left_gatherer.degree(node) + right_gatherer.degree(node);
}


unsigned int Traverser::degree_left(const Kmer& node) const
{
    return left_gatherer.degree(node);
}


unsigned int Traverser::degree_right(const Kmer& node) const
{
    return right_gatherer.degree(node);
}




/******************************************
 * AssemblerTraverser
 ******************************************/

template<bool direction>
AssemblerTraverser<direction>::AssemblerTraverser(const Hashgraph * ht,
                                                  Kmer start_kmer,
                                                  KmerFilterList filters) :
        NodeCursor<direction>(ht, start_kmer, filters)
{
    visited = std::make_shared<SeenSet>();
    AssemblerTraverser<direction>::push_filter(get_visited_filter(visited));
}

template<bool direction>
AssemblerTraverser<direction>::AssemblerTraverser(const Hashgraph * ht,
                                                  Kmer start_kmer,
                                                  KmerFilterList filters,
                                                  std::shared_ptr<SeenSet> visited) :
        NodeCursor<direction>(ht, start_kmer, filters), visited(visited)
{
    AssemblerTraverser<direction>::push_filter(get_visited_filter(visited));
}

template<bool direction>
AssemblerTraverser<direction>::AssemblerTraverser(const AssemblerTraverser<direction>& other) : 
    AssemblerTraverser<direction>(other.graph,
                                  other.cursor,
                                  other.filters,
                                  other.visited)
{
}

template <>
std::string AssemblerTraverser<TRAVERSAL_RIGHT>::join_contigs(std::string& contig_a,
        std::string& contig_b, WordLength offset)
const
{
    return contig_a + contig_b.substr(_ksize - offset);
}

template <>
std::string AssemblerTraverser<TRAVERSAL_LEFT>::join_contigs(std::string& contig_a,
        std::string& contig_b, WordLength offset)
const
{
    return contig_b + contig_a.substr(_ksize - offset);
}

template<bool direction>
char AssemblerTraverser<direction>::next_symbol()
{
    short found = 0;
    char found_base = '\0';
    Kmer neighbor;
    Kmer cursor_next;

    visited->insert(this->cursor);
    for (auto base : alphabets::DNA_SIMPLE) {
        // Get the putative neighbor for this base at the cursor position
        neighbor = NodeCursor<direction>::get_neighbor(this->cursor, base);

        // Now check that the putative neighbor is in the graph and passes the filters
        if (this->graph->get_count(neighbor) &&
                !apply_kmer_filters(neighbor, this->filters)) {

            found++;
            // This naive traverser stops on high degree nodes
            if (found > 1) {
                return '\0';
            }
            found_base = base;
            cursor_next = neighbor;
        }
    }

    if (!found) {
        return '\0';
    } else {
        this->cursor = cursor_next;
        return found_base;
    }
}


template class NodeGatherer<TRAVERSAL_LEFT>;
template class NodeGatherer<TRAVERSAL_RIGHT>;
template class NodeCursor<TRAVERSAL_LEFT>;
template class NodeCursor<TRAVERSAL_RIGHT>;
template class AssemblerTraverser<TRAVERSAL_RIGHT>;
template class AssemblerTraverser<TRAVERSAL_LEFT>;


} // namespace oxli
