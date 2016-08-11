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

template<bool direction>
AssemblerTraverser<direction>::AssemblerTraverser(const Hashtable * ht,
                                 Kmer start_kmer,
                                 KmerFilterList filters) :
    Traverser(ht), filters(filters)
{
    cursor = start_kmer;
}

template<>
Kmer AssemblerTraverser<LEFT>::get_neighbor(Kmer& node,
                                                 const char symbol) {
    return get_left(node, symbol);
}

template<>
Kmer AssemblerTraverser<RIGHT>::get_neighbor(Kmer& node,
                                             const char symbol) {
    return get_right(node, symbol);
}

template<bool direction>
char AssemblerTraverser<direction>::next_symbol()
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
                return '\0';
            }
            base = *symbol_ptr;
            cursor_next = neighbor;
        }
        symbol_ptr++;
    }
    #if DEBUG
    std::cout << "** end next_symbol: " << found << " neighbors. **" << std::endl;
    #endif
    if (!found) {
        return '\0';
    } else {
        cursor = cursor_next;
        return base;
    }
}

template<bool direction>
bool AssemblerTraverser<direction>::set_cursor(Kmer& node)
{
    if(!apply_kmer_filters(node, filters)) {
        cursor = node;
        return true;
    }
    return false;
}


template<bool direction>
Kmer AssemblerTraverser<direction>::get_cursor()
{
    return cursor;
}

template<bool direction>
void AssemblerTraverser<direction>::add_filter(KmerFlter filter)
{
    node_filters.push_back(filter);
}


LinearAssembler::LinearAssembler(const Hashtable * ht) :
    graph(ht), _ksize(ht->ksize())
{

}

KmerFilter LinearAssembler::get_stop_bf_filter(const Hashtable * stop_bf) const
{
    auto filter = [&] (Kmer& n) {
        return stop_bf->get_count(n);
    };
    return filter;
}


// Starting from the given seed k-mer, assemble the maximal linear path in
// both directions.
//
// No guarantees on direction, of course - this may return the reverse
// complement of the input sequence.
std::string LinearAssembler::assemble(const Kmer seed_kmer,
                                      const Hashtable * stop_bf)
    const
{
    std::string right_contig = assemble_right(seed_kmer, stop_bf);
    std::string left_contig = assemble_left(seed_kmer, stop_bf);

    #if DEBUG
    std::cout << "Left: " << left_contig << std::endl;
    std::cout << "Right: " << right_contig << std::endl;
    #endif

    right_contig = right_contig.substr(_ksize);
    return left_contig + right_contig;
}


std::string LinearAssembler::assemble_right(const Kmer seed_kmer,
                                            const Hashtable * stop_bf)
    const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<RIGHT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed(cursor);
}


std::string LinearAssembler::assemble_left(const Kmer seed_kmer,
                                      const Hashtable * stop_bf)
    const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<LEFT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed(cursor);
}


std::string LinearAssembler::_assemble_directed(AssemblerTraverser<LEFT>& cursor)
    const
{
    std::string contig = cursor.get_cursor().get_string_rep(_ksize);
    if (!cursor.get_cursor().is_forward()) {
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

    return contig;
}


std::string LinearAssembler::_assemble_directed(AssemblerTraverser<RIGHT>& cursor)
    const
{
    std::string contig = cursor.get_cursor().get_string_rep(_ksize);
    char next_base;

    #if DEBUG
    std::cout << "## assemble_right\nContig: " << contig << std::endl;
    #endif

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
    }

    return contig;
}



/********************************
 * Labeled Assembly
 ********************************/

LabeledLinearAssembler::LabeledLinearAssembler(const LabelHash * labels) :
    graph(labels->graph), lh(lh), _ksize(graph->ksize())
{

}


// Starting from the given seed k-mer, assemble all maximal linear paths in
// both directions, using labels to skip over tricky bits.

StringVector LabeledLinearAssembler::assemble(const Kmer seed_kmer,
                                              const Hashtable * stop_bf=0)
    const
{
    KmerFilterList node_filters;
    if (stop_bf) {
        node_filters.push_back(linear_asm->get_stop_bf_filter(stop_bf);
    }

    SeenSet visited;
    KmerFilter visited_filter = [&] (Kmer& n) {
        if(set_contains(visited, n)) {
            return true;
        } else {
            visited.insert(n);
            return false;
        }
    };
    node_filters.push_back(visited_filter);

#if DEBUG
    std::cout << "assemble right: " << start_kmer << std::endl;
#endif

    StringVector fwd_paths;
    AssemblerTraverser<RIGHT> rcursor(graph, seed_kmer, node_filters);
    _assemble_directed(rcursor, fwd_paths);
/*
#if DEBUG
    std::cout << "assemble left: " << start_kmer << std::endl;
#endif

    start_kmer = _revcomp(start_kmer);
    std::vector<std::string> rev_paths;
    _assemble_labeled_right(start_kmer.c_str(), rev_paths, visited);
    visited.clear();
#if DEBUG
    std::cout << "join right and left contigs: " << rev_paths.size() << std::endl;
#endif
    std::vector<std::string> paths;
    for (unsigned int i = 0; i < rev_paths.size(); i++) {
        for (unsigned int j = 0; j < fwd_paths.size(); j++) {
            std::string left = rev_paths[i];
            left = left.substr(graph->_ksize);
            std::string contig = _revcomp(left) + fwd_paths[j];
            paths.push_back(contig);
        }
    }
*/
    return fwd_paths;
}

void LabeledLinearAssembler::_assemble_directed(AssemblerTraverser<RIGHT>& start_cursor,
                                                StringVector& paths)
    const
{

    std::string root_contig = linear_asm->_assemble_directed(start_cursor);
    Kmer end_kmer = start_cursor.get_cursor();

    if (right_degree(end_kmer) > 1) {               // hit a HDN
#if DEBUG
        std::cout << "Contig thus far: " << root_contig << std::endl;
        std::cout << "HDN: " << end_kmer.repr(_ksize) << "\n";
#endif // DEBUG

        LabelSet labels;
        lh->get_tag_labels(end_kmer, labels);

        /* For each label, we try to find spanning paths. We create
         * a new cursor starting at the end k-mer, with the existing node
         * filters; then we give that to the label spanning function.
         *
         * NOTE: This implies that there may be some non-deterministic
         * behavior: the ordering of the labels could vary, and decides
         * the recursion termination.
         */
        if(labels.size() == 0) {
            // if no labels are found there's nothing to be done, return
            paths.push_back(root_contig);
            return;
        } else {
            for (auto label : labels) {
                // Copy the current cursor at end_cursor for the spanning function
                AssemblerTraverser<RIGHT> span_cursor(start_cursor);
                std::string spanning_contig = assemble_across_label(span_cursor, label);
                std::string extended_contig = root_contig

                // Creata a new cursor with the start cursor's filters
                // but set it to the end of the extended contig
                AssemblerTraverser<RIGHT> continue_cursor(start_cursor);
                continue_cursor.set_cursor(span_cursor.get_cursor());

                // Recurse and gather paths
                StringVector continue_contigs;
                _assemble_directed(continue_cursor, continue_contigs);

                if (continue_contigs.size() == 0) {
                    paths.push_back(extended_contig);
                } else {
                    for (auto continue_contig : continue_contigs) {
                        std::string full_contig = root_contig + spanning_contig.substr(_ksize);
                        full_contig += continue_contig.substr(_ksize);
                        paths.push_back(full_contig);
                    }
                }
            } //end for
        }
    } else {
        paths.push_back(contig);
    }
}

//TODO make member function to create this filter
std::string LabeledLinearAssembler::_assemble_across_labels(AssemblerTraverser& start_cursor,
                                                            const Label label)
    const
{
    /* We'll now assemble, following the given label, as far as we can.
     * We add an extra filter to the list: now, if we find no labels, we
     * continue assembling; if we find labels and ours is included, we
     * continue; and if we find labels and ours is not included, we stop.
     * This cursor should also have the filters for visited k-mers and
     * the stop bloom filter already.
     */
    KmerFilter label_filter = [&] (Kmer& node) {
        LabelSet ls;
        lh->get_tag_labels(node, ls);
        if (ls.size() == 0) {
            return false;
        } else {
            return !set_contains(ls, label);
        }
    };
    start_cursor.add_filter(label_filter);

    return linear_asm._assemble_directed(start_cursor);
}
