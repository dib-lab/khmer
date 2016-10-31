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

#include "assembler.hh"

#include <algorithm>
#include <iostream>

using namespace std;

namespace khmer
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
    std::string right_contig = assemble_right(seed_kmer, stop_bf);
    std::string left_contig = assemble_left(seed_kmer, stop_bf);

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

    AssemblerTraverser<RIGHT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<RIGHT>(cursor);
}


std::string LinearAssembler::assemble_left(const Kmer seed_kmer,
        const Hashgraph * stop_bf)
const
{
    std::list<KmerFilter> node_filters;
    if (stop_bf) {
        node_filters.push_back(get_stop_bf_filter(stop_bf));
    }

    AssemblerTraverser<LEFT> cursor(graph, seed_kmer, node_filters);
    return _assemble_directed<LEFT>(cursor);
}

template <>
std::string LinearAssembler::_assemble_directed<LEFT>(AssemblerTraverser<LEFT>&
        cursor)
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
std::string LinearAssembler::_assemble_directed<RIGHT>
(AssemblerTraverser<RIGHT>& cursor)
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

    SeenSet * visited = new SeenSet();

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled RIGHT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector right_paths;
    NonLoopingAT<RIGHT> rcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<RIGHT>(rcursor, right_paths);

#if DEBUG_ASSEMBLY
    std::cout << "Assemble Labeled LEFT: " << seed_kmer.repr(_ksize) << std::endl;
#endif
    StringVector left_paths;
    NonLoopingAT<LEFT> lcursor(graph, seed_kmer, node_filters, visited);
    _assemble_directed<LEFT>(lcursor, left_paths);

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
void SimpleLabeledAssembler::_assemble_directed(NonLoopingAT<direction>&
        start_cursor,
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
    std::vector< NonLoopingAT<direction> > cursors;

    segments.push_back(root_contig);
    cursors.push_back(start_cursor);

    while(segments.size() != 0) {

        std::string segment = segments.back();
        NonLoopingAT<direction> cursor = cursors.back();
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
                    NonLoopingAT<direction> branch_cursor(cursor);
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



/*
template <bool direction>
void SimpleLabeledAssembler::_assemble_directed(NonLoopingAT<direction>&
        start_cursor,
        StringVector& paths)
const
{
#if DEBUG_ASSEMBLY
    std::cout << "## assemble_labeled_directed_" << direction << " [start] at " <<
        start_cursor.cursor.repr(_ksize) << std::endl;
#endif
    std::string root_contig = linear_asm->_assemble_directed<direction>
                              (start_cursor);

    if (start_cursor.cursor_degree() > 1) {               // hit a HDN
#if DEBUG_ASSEMBLY
        std::cout << "Root contig: " << root_contig << std::endl;
        std::cout << "HDN: " << start_cursor.cursor.repr(_ksize) << "\n";
        std::cout << start_cursor.cursor_degree() << std::endl;
#endif // DEBUG_ASSEMBLY

        LabelSet labels;
        lh->get_tag_labels(start_cursor.cursor, labels);


        if(labels.size() == 0) {
            // if no labels are found there's nothing to be done, return

            paths.push_back(root_contig);
            return;
        } else {
#if DEBUG_ASSEMBLY
            std::cout << "Found " << labels.size() << " labels" << std::endl;
#endif



            start_cursor.push_filter(get_simple_label_intersect_filter(labels, lh));
            KmerQueue branch_starts;
            start_cursor.neighbors(branch_starts);
            start_cursor.pop_filter();
#if DEBUG_ASSEMBLY
                std::cout << branch_starts.size() << " neighbors found" << std::endl;
#endif
            if (branch_starts.size() == 0) {
                paths.push_back(root_contig);
                return;
            }

            StringVector branch_contigs;
            while(!branch_starts.empty()) { // TODO: change from queue
                NonLoopingAT<direction> branch_cursor(start_cursor);
                branch_cursor.cursor = branch_starts.front();
                branch_starts.pop();

                _assemble_directed<direction>(branch_cursor, branch_contigs);
            }

            for (auto branch_contig : branch_contigs) {
                std::string full_contig = start_cursor.join_contigs(root_contig,
                                                                   branch_contig,
                                                                   1);
                paths.push_back(full_contig);
            }
        }
    } else {
        paths.push_back(root_contig);
    }
}

*/

}
