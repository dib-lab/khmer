/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2015, The Regents of the University of California.

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
#include "hashtable.hh"
#include "traversal.hh"

using namespace khmer;
using namespace std;

Traverser::Traverser(const Hashtable * ht) :
    KmerFactory(ht->ksize()), graph(ht)
{
    bitmask = 0;
    for (unsigned int i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    rc_left_shift = _ksize * 2 - 2;
}

Kmer Traverser::get_left(Kmer& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = ((node.kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift);
    kmer_r = (((node.kmer_r) << 2) & bitmask) | (twobit_comp(ch));
    return build_kmer(kmer_f, kmer_r);
}


Kmer Traverser::get_right(Kmer& node, const char ch)
{
    HashIntoType kmer_f, kmer_r;
    kmer_f = (((node.kmer_f) << 2) & bitmask) | (twobit_repr(ch));
    kmer_r = ((node.kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift);
    return build_kmer(kmer_f, kmer_r);
}

unsigned int Traverser::traverse_left(Kmer& node,
                                      KmerQueue & node_q,
                                      std::function<bool (Kmer&)> filter)
{
    unsigned int found = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node) && (!filter || filter(prev_node))) {
            node_q.push(prev_node);
            ++found;
        }
        ++base;
    }

    return found;
}

unsigned int Traverser::traverse_right(Kmer& node,
                                       KmerQueue & node_q,
                                       std::function<bool (Kmer&)> filter)
{
    unsigned int found = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node) && (!filter || filter(next_node))) {
            node_q.push(next_node);
            ++found;
        }
        ++base;
    }

    return found;
}

unsigned int Traverser::degree_left(Kmer& node)
{
    unsigned int degree = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer prev_node = get_left(node, *base);
        if (graph->get_count(prev_node)) {
            ++degree;
        }
        ++base;
    }

    return degree;
}

unsigned int Traverser::degree_right(Kmer& node)
{
    unsigned int degree = 0;

    char bases[] = "ACGT";
    char * base = bases;
    while(*base != '\0') {
        Kmer next_node = get_right(node, *base);
        if (graph->get_count(next_node)) {
            ++degree;
        }
        ++base;
    }

    return degree;
}

unsigned int Traverser::degree(Kmer& node)
{
    return degree_right(node) + degree_left(node);
}
