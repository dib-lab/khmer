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
#include "hashtable.hh"
#include "traversal.hh"
#include "assembler.hh"
#include "symbols.hh"

#include <algorithm>
#include <iostream>

using namespace khmer;
using namespace std;

AssemblerTraverser::AssemblerTraverser(const Hashtable * ht,
                                 Kmer start_kmer,
                                 KmerFilterList filters,
                                 bool traverse_right) :
    Traverser(ht), filters(filters), traverse_right(traverse_right)
{
    cursor = start_kmer;
    contig_kmers.push_back(cursor);
}

void AssemblerTraverser::gather_linear_path()
{
    while (1) {
        char * base = alphabets::DNA_SIMPLE;
        short found = 0;
        Kmer neighbor;

        while(*base != '\0') {

            if(traverse_right) { // NOTE: hoping this gets optimized out because const
                neighbor = get_right(cursor, *base);
            } else {
                neighbor = get_left(cursor, *base);
            }
            std::cout << "Try: " << (char)*base << " (" << neighbor << ")" << std::endl;
            if (graph->get_count(neighbor) &&
                !apply_kmer_filters(neighbor, filters)) {
                std::cout << "Found " << (char)*base << std::endl;
                found++;
                if (found > 1) {
                    break;
                }
            }
            base++;
        }
        std::cout << "Found: " << found << std::endl;
        if (found != 1) {
            break;
        } else {
            contig_kmers.push_back(neighbor);
            cursor = neighbor;
        }
    }
}

unsigned int AssemblerTraverser::get_path_length()
    const
{
    return _ksize + (contig_kmers.size() - 1);
}

std::string AssemblerTraverser::build_contig()
    const
{

    if (traverse_right) {
        auto it = contig_kmers.begin();
        std::string contig = (*it).get_string_rep(_ksize);
        ++it;

        for (; it!=contig_kmers.end(); it++) {
            //contig += revtwobit_repr((*it) & 3); // just get the last base
            contig += (*it).get_last_base();
        }
        std::cout << "Right (" << contig.length() << "):" << contig << std::endl;
        std::cout << "Path Length: " << get_path_length() << std::endl;
        return contig;
    } else {
        auto it = contig_kmers.rbegin();
        std::string contig = (*it).get_string_rep(_ksize);
        ++it;

        for (; it!=contig_kmers.rend(); it++) {
            //contig += revtwobit_repr((*it) & 3); // just get the last base
            contig += (*it).get_last_base();

        }
        std::cout << "Left (" << contig.length() << "):" << contig << std::endl;
        std::cout << "Left Length: " << get_path_length() << std::endl;
        return contig;
    }
}

std::string AssemblerTraverser::assemble()
{
    gather_linear_path();
    return build_contig();
}

Assembler::Assembler(const Hashtable * ht) :
    Traverser(ht)
{

}


// Starting from the given seed k-mer, assemble the maximal linear path in
// both directions.
//
// No guarantees on direction, of course - this may return the reverse
// complement of the input sequence.
//
// Note: as written, will ignore branches to the left and continue
// past them; this probably needs to be fixed.  For now, this means
// that assembling from two different directions may yield different
// results.
std::string Assembler::assemble_linear_path(const Kmer seed_kmer,
                                            const Hashtable * stop_bf)
    const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        auto stop_bf_filter = [&] (Kmer& n) {
            return stop_bf->get_count(n);
        };
        node_filters.push_back(stop_bf_filter);
    }

    AssemblerTraverser right_traverser(graph, seed_kmer, node_filters);
    AssemblerTraverser left_traverser(graph, seed_kmer, node_filters, 0);

    //std::string start_kmer = seed_kmer.get_string_rep(_ksize);

    //std::string right = _assemble_right(start_kmer, node_filters);
    //std::string left = _assemble_left(start_kmer, node_filters);

    std::string right = right_traverser.assemble();
    std::string left = left_traverser.assemble();

    right = right.substr(_ksize);
    return left + right;
}


std::string Assembler::_assemble_left(const std::string start_kmer,
                                      std::list<KmerFilter>& node_filters)
    const
{
    std::string contig = _assemble_right(_revcomp(start_kmer), node_filters);
    return _revcomp(contig);
}


std::string Assembler::_assemble_right(const std::string start_kmer,
                                       std::list<KmerFilter>& node_filters)
    const
{
    std::string kmer = start_kmer;
    std::string contig = kmer;

    while (1) {
        char * base = alphabets::DNA_SIMPLE;
        bool found = false;
        char found_base;
        bool found2 = false;

        while(*base != 0) {
            std::string try_kmer = kmer.substr(1) + (char) *base;
            Kmer try_hashed = build_kmer(try_kmer);

            // a hit!
            if (graph->get_count(try_hashed) &&
                !apply_kmer_filters(try_hashed, node_filters)) {

                if (found) {
                    found2 = true;
                    break;
                }
                found_base = (char) *base;
                found = true;
            }
            base++;
        }
        if (!found || found2) {
            break;
        } else {
            contig += found_base;
            kmer = kmer.substr(1) + found_base;
            found = true;
        }
    }
    return contig;
}
