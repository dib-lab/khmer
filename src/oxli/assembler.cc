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

#include "oxli/assembler.hh"

#include <algorithm>
#include <iostream>

using namespace std;

namespace oxli
{

/********************************
 * Simple Linear Assembly
 ********************************/

LinearAssembler::LinearAssembler(const Hashgraph * ht) :
    graph(ht), _ksize(ht->ksize())
{

}

// Starting from the given seed k-mer, assemble the maximal linear path in
// both directions.
std::string LinearAssembler::assemble(const Kmer seed_kmer,
                                      const Hashgraph * stop_bf)
const
{
    if (graph->get_count(seed_kmer) == 0) {
        // If the seed k-mer is not in the de Bruijn graph, stop trying to make
        // something happen. It's not going to happen!
        return "";
    }

    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    std::shared_ptr<SeenSet> visited = std::make_shared<SeenSet>();
    AssemblerTraverser<TRAVERSAL_RIGHT> rcursor(graph, seed_kmer, node_filters, visited);
    AssemblerTraverser<TRAVERSAL_LEFT> lcursor(graph, seed_kmer, node_filters, visited);

    std::string right_contig = _assemble_directed<TRAVERSAL_RIGHT>(rcursor);
    std::string left_contig = _assemble_directed<TRAVERSAL_LEFT>(lcursor);

#if DEBUG_ASSEMBLY
    std::cout << "Left: " << left_contig << std::endl;
    std::cout << "Right: " << right_contig << std::endl;
#endif

    right_contig = right_contig.substr(_ksize);
    return left_contig + right_contig;
}


std::string LinearAssembler::assemble_right(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<TRAVERSAL_RIGHT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<TRAVERSAL_RIGHT>(cursor);
}


std::string LinearAssembler::assemble_left(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<TRAVERSAL_LEFT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<TRAVERSAL_LEFT>(cursor);
}

template <>
std::string LinearAssembler::
_assemble_directed<TRAVERSAL_LEFT>(AssemblerTraverser<TRAVERSAL_LEFT>& cursor)
const
{
    std::string contig = cursor.cursor.get_string_rep(_ksize);
    if (!cursor.cursor.is_forward()) {
        contig = _revcomp(contig);
    }

#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_left[start] at " << contig << std::endl;
#endif

    reverse(contig.begin(), contig.end());
    char next_base;
    unsigned int found = 0;

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
        found++;
    }

    reverse(contig.begin(), contig.end());
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_left[end] found " << found << std::endl;
#endif

    return contig;
}

template<>
std::string LinearAssembler::
_assemble_directed<TRAVERSAL_RIGHT>(AssemblerTraverser<TRAVERSAL_RIGHT>& cursor)
const
{
    std::string contig = cursor.cursor.get_string_rep(_ksize);
    if (!cursor.cursor.is_forward()) {
        contig = _revcomp(contig);
    }
    char next_base;
    unsigned int found = 0;

#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_right[start] at " << contig << std::endl;
#endif

    while ((next_base = cursor.next_symbol()) != '\0') {
        contig += next_base;
        found++;
    }
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_linear_right[end] found " << found << std::endl;
#endif
    return contig;
}


/********************************
 * Labeled Assembly
 ********************************/

SimpleLabeledAssembler::SimpleLabeledAssembler(const LabelHash * lh) :
    graph(lh->graph), lh(lh), _ksize(lh->graph->ksize())
{
    linear_asm = new LinearAssembler(graph);
}

SimpleLabeledAssembler::~SimpleLabeledAssembler()
{
    delete this->linear_asm;
}


// Starting from the given seed k-mer, assemble all maximal linear paths in
// both directions, using labels to skip over tricky bits.
StringVector SimpleLabeledAssembler::assemble(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled: " << seed_kmer.repr(_ksize) << std::endl;
#endif

    KmerFilterList node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    std::shared_ptr<SeenSet> visited = std::make_shared<SeenSet>();

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled RIGHT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector right_paths;
    AssemblerTraverser<TRAVERSAL_RIGHT> rcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<TRAVERSAL_RIGHT>(rcursor, right_paths);

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled LEFT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector left_paths;
    AssemblerTraverser<TRAVERSAL_LEFT> lcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<TRAVERSAL_LEFT>(lcursor, left_paths);

    StringVector paths;
    for (unsigned int i = 0; i < left_paths.size(); i++) {
        for (unsigned int j = 0; j < right_paths.size(); j++) {
            std::string right = right_paths[j];
            right = right.substr(_ksize);
            std::string contig = left_paths[i] + right;
            paths.push_back(contig);
        }
    }

    visited->clear();
    return paths;
}

template <bool direction>
void SimpleLabeledAssembler::
_assemble_directed(AssemblerTraverser<direction>& start_cursor,
        		   StringVector& paths)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_labeled_directed_" << direction << " [start] at " <<
              start_cursor.cursor.repr(_ksize) << std::endl;
#endif

    // prime the traversal with the first linear segment
    std::string root_contig = linear_asm->_assemble_directed<direction>
                              (start_cursor);
#if DEBUG_ASSEMBLY
    std::cout << "Primed: " << root_contig << std::endl;
    std::cout << "Cursor: " << start_cursor.cursor.repr(_ksize) << std::endl;
#endif
    StringVector segments;
    std::vector< AssemblerTraverser<direction> > cursors;

    segments.push_back(root_contig);
    cursors.push_back(start_cursor);

    while(segments.size() != 0) {

        std::string segment = segments.back();
        AssemblerTraverser<direction> cursor = cursors.back();
#if DEBUG_ASSEMBLY
        std::cout << "Pop: " << segments.size() << " segments on stack." << std::endl;
        std::cout << "Segment: " << segment << std::endl;
        std::cout << "Cursor: " << cursor.cursor.repr(_ksize) << std::endl;
        std::cout << "n_filters: " << cursor.n_filters() << std::endl;
#endif
        segments.pop_back();
        cursors.pop_back();

        // check if the cursor has hit a HDN or reached a dead end
        if (cursor.cursor_degree() > 1) {

            LabelSet labels;
            lh->get_tag_labels(cursor.cursor, labels);

            if(labels.size() == 0) {
                // if no labels are found we can do nothing; gather this contig
#if DEBUG_ASSEMBLY
                std::cout << "no-label dead-end" << std::endl;
#endif
                paths.push_back(segment);
                continue;
            } else {
                // if there are labels, try to hop the HDN.
                // first, get a label filter
                cursor.push_filter(get_simple_label_intersect_filter(labels, lh));
                KmerQueue branch_starts;
                // now get neighbors that pass the filter
                cursor.neighbors(branch_starts);
                // remove the filter
                cursor.pop_filter();

                // no neighbors found; done with this path
                if (branch_starts.empty()) {
#if DEBUG_ASSEMBLY
                    std::cout << "no-neighbors dead-end" << std::endl;
#endif
                    paths.push_back(segment);
                    continue;
                }

                // found some neighbors; extend them
                while(!branch_starts.empty()) {
                    // spin off a cursor for the new branch
                    AssemblerTraverser<direction> branch_cursor(cursor);
                    branch_cursor.cursor = branch_starts.front();
                    branch_starts.pop();

#if DEBUG_ASSEMBLY
                    std::cout << "Branch cursor: " << branch_cursor.cursor.repr(
                                  _ksize) << std::endl;
#endif

                    // assemble linearly as far as possible
                    std::string branch = linear_asm->_assemble_directed<direction>(branch_cursor);
                    // create a new segment with the branch
                    std::string new_segment = branch_cursor.join_contigs(segment, branch, 1);
#if DEBUG_ASSEMBLY
                    std::cout << "Push segment: " << new_segment << std::endl;
                    std::cout << "Push cursor: " << branch_cursor.cursor.repr(_ksize) << std::endl;
#endif
                    segments.push_back(new_segment);
                    cursors.push_back(branch_cursor);
                }
            }
        } else {
            // this segment is a dead-end; keep the contig
#if DEBUG_ASSEMBLY
            std::cout << "degree-1 dead-end" << std::endl;
#endif
            paths.push_back(segment);
            continue;
        }
    }
}

/***************************************
 * Junction-counting assembler
 ***************************************/

JunctionCountAssembler::JunctionCountAssembler(Hashgraph * ht) :
    graph(ht), _ksize(ht->ksize()), traverser(ht), linear_asm(ht)
{
    std::vector<uint64_t> table_sizes = graph->get_tablesizes();
    junctions = new Countgraph(_ksize, table_sizes);
}


JunctionCountAssembler::~JunctionCountAssembler()
{
    delete this->junctions;
}

uint16_t JunctionCountAssembler::consume(std::string sequence)
{
    // first we need to put the sequence in the graph
    graph->consume_string(sequence);
    // now we find its high degree nodes and count the
    // branch junctions
    KmerIterator kmers(sequence.c_str(), _ksize);
    Kmer kmer = kmers.next();
    if (kmers.done()) {
        return 0;
    }
    Kmer next_kmer = kmers.next();
    if (kmers.done()) {
        return 0;
    }
    uint16_t d = this->traverser.degree(kmer);
    uint16_t next_d = this->traverser.degree(next_kmer);
    uint16_t n_junctions = 0;

    while(!kmers.done()) {
        if (d > 2 || next_d > 2) {
            count_junction(kmer, next_kmer);
            n_junctions++;
#if DEBUG_ASSEMBLY
            std::cout << "Junction: " << kmer.repr(_ksize) << ", " << next_kmer.repr(
                          _ksize) << std::endl;
            std::cout << "Junction Count: " << get_junction_count(kmer,
                      next_kmer) << std::endl;
#endif
        }
        kmer = next_kmer;
        d = next_d;
        next_kmer = kmers.next();
        next_d = this->traverser.degree(next_kmer);
    }

    return n_junctions / 2;
}

void JunctionCountAssembler::count_junction(Kmer kmer_a, Kmer kmer_b)
{
    junctions->count(kmer_a.kmer_u ^ kmer_b.kmer_u);
}

BoundedCounterType JunctionCountAssembler::get_junction_count(Kmer kmer_a,
        Kmer kmer_b)
const
{
    return junctions->get_count(kmer_a.kmer_u ^ kmer_b.kmer_u);
}

// Starting from the given seed k-mer, assemble all maximal linear paths in
// both directions, using junction counts to skip over tricky bits.
StringVector JunctionCountAssembler::assemble(const Kmer seed_kmer,
                                              const Hashtable * stop_bf)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "Assemble Junctions: " << seed_kmer.repr(_ksize) << std::endl;
#endif

    KmerFilterList node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    std::shared_ptr<SeenSet> visited = std::make_shared<SeenSet>();

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Junctions RIGHT: " << seed_kmer.repr(
                  _ksize) << std::endl;
#endif
    StringVector right_paths;
    AssemblerTraverser<TRAVERSAL_RIGHT> rcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<TRAVERSAL_RIGHT>(rcursor, right_paths);

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Junctions LEFT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector left_paths;
    AssemblerTraverser<TRAVERSAL_LEFT> lcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<TRAVERSAL_LEFT>(lcursor, left_paths);

    StringVector paths;
    for (unsigned int i = 0; i < left_paths.size(); i++) {
        for (unsigned int j = 0; j < right_paths.size(); j++) {
            std::string right = right_paths[j];
            right = right.substr(_ksize);
            std::string contig = left_paths[i] + right;
            paths.push_back(contig);
        }
    }

    visited->clear();
    return paths;
}

template <bool direction>
void JunctionCountAssembler::
_assemble_directed(AssemblerTraverser<direction>& start_cursor,
        		   StringVector& paths)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_junctions_directed_" << direction << " [start] at " <<
              start_cursor.cursor.repr(_ksize) << std::endl;
#endif

    // prime the traversal with the first linear segment
    std::string root_contig = linear_asm._assemble_directed<direction>
                              (start_cursor);
#if DEBUG_ASSEMBLY
    std::cout << "Primed: " << root_contig << std::endl;
    std::cout << "Cursor: " << start_cursor.cursor.repr(_ksize) << std::endl;
#endif
    StringVector segments;
    std::vector< AssemblerTraverser<direction> > cursors;

    segments.push_back(root_contig);
    cursors.push_back(start_cursor);

    while(segments.size() != 0) {

        std::string segment = segments.back();
        AssemblerTraverser<direction> cursor = cursors.back();
#if DEBUG_ASSEMBLY
        std::cout << "Pop: " << segments.size() << " segments on stack." << std::endl;
        std::cout << "Segment: " << segment << std::endl;
        std::cout << "Cursor: " << cursor.cursor.repr(_ksize) << std::endl;
        std::cout << "n_filters: " << cursor.n_filters() << std::endl;
#endif
        segments.pop_back();
        cursors.pop_back();

        // check if the cursor has hit a HDN or reached a dead end
        if (cursor.cursor_degree() > 1) {

            cursor.push_filter(get_junction_count_filter(cursor.cursor, this->junctions));
            KmerQueue branch_starts;
            // now get neighbors that pass the filter
            cursor.neighbors(branch_starts);
            // remove the filter
            cursor.pop_filter();

            // no neighbors found; done with this path
            if (branch_starts.empty()) {
                paths.push_back(segment);
                continue;
            }

            // found some neighbors; extend them
            while(!branch_starts.empty()) {
                // spin off a cursor for the new branch
                AssemblerTraverser<direction> branch_cursor(cursor);

                branch_cursor.cursor = branch_starts.front();
                branch_starts.pop();

                // assemble linearly as far as possible
                std::string branch = linear_asm._assemble_directed<direction>(branch_cursor);
                // create a new segment with the branch
                std::string new_segment = branch_cursor.join_contigs(segment, branch, 1);

                segments.push_back(new_segment);
                cursors.push_back(branch_cursor);
            }
        } else {
            // this segment is a dead-end; keep the contig
            paths.push_back(segment);
            continue;
        }
    }
}

}
