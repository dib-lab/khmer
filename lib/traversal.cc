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
#include "khmer.hh"
#include "hashtable.hh"
#include "traversal.hh"
#include "alphabets.hh"
#include "kmer_hash.hh"

#define DEBUG 1

using namespace std;

namespace khmer
{

/******************************************
 * NodeGatherer
 ******************************************/

template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashtable * ht,
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
NodeGatherer<direction>::NodeGatherer(const Hashtable * ht) : 
    NodeGatherer(ht, KmerFilterList())
{
}


template <bool direction>
NodeGatherer<direction>::NodeGatherer(const Hashtable * ht, 
                                      KmerFilter filter) : 
    NodeGatherer(ht, KmerFilterList())
{
    filters.push_back(filter);
}


template<>
Kmer NodeGatherer<LEFT>::get_neighbor(const Kmer& node, const char ch)
const
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return build_kmer(kmer_f, kmer_r);
}


template<>
Kmer NodeGatherer<RIGHT>::get_neighbor(const Kmer& node,
                                       const char ch)
const
{
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
        Kmer neighbor = get_neighbor(node, base);
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
NodeCursor<direction>::NodeCursor(const Hashtable * ht,
                                Kmer start_kmer,
                                KmerFilterList filters) :
    NodeGatherer<direction>(ht, filters)
{
    cursor = start_kmer;
}


template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashtable * ht,
                                Kmer start_kmer) :
    NodeCursor<direction>(ht, start_kmer, KmerFilterList())
{
}


template<bool direction>
NodeCursor<direction>::NodeCursor(const Hashtable * ht,
                                  Kmer start_kmer,
                                  KmerFilter filter) :
    NodeCursor<direction>(ht, start_kmer)
{
    push_filter(filter);
}


/******************************************
 * Traverser
 ******************************************/

Traverser::Traverser(const Hashtable * ht,
                     KmerFilterList filters) :
    KmerFactory(ht->ksize()), 
    graph(ht), 
    left_gatherer(ht, filters),
    right_gatherer(ht, filters)
{
}

Traverser::Traverser(const Hashtable * ht, 
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
unsigned int AssemblerTraverser<direction>::cursor_degree()
const
{
    return this->degree(this->cursor);
}


template <>
std::string AssemblerTraverser<RIGHT>::join_contigs(std::string& contig_a,
        std::string& contig_b, WordLength offset)
const
{
    return contig_a + contig_b.substr(_ksize - offset);
}

template <>
std::string AssemblerTraverser<LEFT>::join_contigs(std::string& contig_a,
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

    for (auto base : alphabets::DNA_SIMPLE) {
        neighbor = NodeCursor<direction>::get_neighbor(this->cursor, base);

        if (this->graph->get_count(neighbor) &&
                !apply_kmer_filters(neighbor, this->filters)) {

            found++;
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

/******************************************
 * NonLoopingAT
 ******************************************/

template<bool direction>
NonLoopingAT<direction>::NonLoopingAT(const Hashtable * ht,
                                      Kmer start_kmer,
                                      KmerFilterList filters,
                                      SeenSet * visited) :
    AssemblerTraverser<direction>(ht, start_kmer, filters), visited(visited)
{
    AssemblerTraverser<direction>::push_filter(get_visited_filter(visited));
}

template<bool direction>
char NonLoopingAT<direction>::next_symbol()
{
    visited->insert(this->cursor);
    return AssemblerTraverser<direction>::next_symbol();
}

template class NodeGatherer<LEFT>;
template class NodeGatherer<RIGHT>;
template class NodeCursor<LEFT>;
template class NodeCursor<RIGHT>;
template class AssemblerTraverser<RIGHT>;
template class AssemblerTraverser<LEFT>;
template class NonLoopingAT<RIGHT>;
template class NonLoopingAT<LEFT>;


} // namespace khmer
