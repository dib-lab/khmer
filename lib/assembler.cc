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

#define DEBUG 0

using namespace khmer;
using namespace std;

<template bool direction>
AssemblerTraverser::AssemblerTraverser(const Hashtable * ht,
                                 Kmer start_kmer,
                                 KmerFilterList filters,
                                 bool direction) :
    Traverser(ht), filters(filters), direction(direction)
{
    cursor = start_kmer;
    if(direction == ASSEMBLE_LEFT) {
        redirector = &AssemblerTraverser::get_left;
    }
}

<template bool direction>
Kmer AssemblerTraverser::get_neighbor(Kmer& node, const char symbol) {
    return redirector(this, node, symbol);
}

<template bool direction>
char AssemblerTraverser::next_symbol()
{
    char * symbol_ptr = alphabets::DNA_SIMPLE;
    char base;
    short found = 0;
    Kmer neighbor;
    Kmer cursor_next;

    #if DEBUG
    std::cout << "** start next_symbol (start: " << cursor.repr(_ksize) <<  ") **" << std::endl;
    #endif

    while(*symbol_ptr != '\0') {
        neighbor = get_neighbor(cursor, *symbol_ptr);

        #if DEBUG
            std::cout << "Try: " << *symbol_ptr << " " << neighbor.repr(_ksize)
                << direction << " Count: " << graph->get_count(neighbor) << std::endl;
        #endif

        if (graph->get_count(neighbor) &&
            !apply_kmer_filters(neighbor, filters)) {

            found++;
            if (found > 1) {
                break;
            }
            base = *symbol_ptr;
            cursor_next = neighbor;
        }
        symbol_ptr++;
    }
    #if DEBUG
    std::cout << "** end next_symbol: " << found << " neighbors. **" << std::endl;
    #endif
    if (found != 1) {
        return '\0';
    } else {
        cursor = cursor_next;
        return base;
    }
}

<template bool direction>
bool AssemblerTraverser::set_cursor(Kmer& node)
{
    if(!apply_kmer_filters(node, filters)) {
        cursor = node;
        return true;
    }
    return false;
}

<template bool direction>
Kmer AssemblerTraverser::get_cursor()
{
    return cursor;
}


LinearAssembler::LinearAssembler(const Hashtable * ht) :
    graph(ht), _ksize(ht->ksize())
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
std::string LinearAssembler::assemble(const Kmer seed_kmer,
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

    std::string right_contig;
    AssemblerTraverser<RIGHT> right_cursor(graph, start_kmer, node_filters);
    assemble_right(seed_kmer, right_contig, node_filters);

    std::string left_contig;
    AssemblerTraverser<LEFT> left_cursor(graph, start_kmer, node_filters);
    assemble_left(seed_kmer, left_contig, left_cursor);

    #if DEBUG
    std::cout << "Left: " << left_contig << std::endl;
    std::cout << "Right: " << right_contig << std::endl;
    #endif

    right_contig = right_contig.substr(_ksize);
    return left_contig + right_contig;
}


Kmer LinearAssembler::assemble_left(std::string& contig,
                                    AssemblerTraverser<RIGHT>& cursor)
    const
{
    contig = start_kmer.get_string_rep(_ksize);
    if (!start_kmer.is_forward()) {
        contig = _revcomp(contig);
    }

    #if DEBUG
    std::cout << "## assemble_left\nStart Contig: " << contig << std::endl;
    #endif

    reverse(contig.begin(), contig.end());
    char next_base;

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
    }

    reverse(contig.begin(), contig.end());

    return cursor.get_cursor();
}


Kmer LinearAssembler::assemble_right(std::string& contig,
                                     AssemblerTraverser<LEFT>& cursor)
    const
{
    contig = start_kmer.get_string_rep(_ksize);
    char next_base;

    #if DEBUG
    std::cout << "## assemble_right\nContig: " << contig << std::endl;
    #endif

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
    }
    return contig;
}
