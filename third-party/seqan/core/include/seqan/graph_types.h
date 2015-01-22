// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_TYPES_H
#define SEQAN_HEADER_GRAPH_TYPES_H

// External / STL
#include <deque>


// Seqan
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

// Basic graph stuff
#include <seqan/graph_types/graph_base.h>
#include <seqan/graph_types/graph_idmanager.h>	// Id manager
#include <seqan/graph_types/graph_edgestump.h>	// EdgeStumps
#include <seqan/graph_types/graph_interface.h>	// Graph metafunctions

// Graph types
#include <seqan/graph_types/graph_impl_directed.h>		// Directed Graph
#include <seqan/graph_types/graph_impl_undirected.h>	// Undirected graph
#include <seqan/graph_types/graph_impl_automaton.h>		// Automaton
#include <seqan/graph_types/graph_impl_wordgraph.h>		// Specialized automaton: Word graph
#include <seqan/graph_types/graph_impl_tree.h>			// Tree
#include <seqan/graph_types/graph_impl_fragment.h>		// Fragment
#include <seqan/graph_types/graph_impl_hmm.h>			// HMM

// Graph iterators
#include <seqan/graph_types/graph_iterator.h>
#include <seqan/graph_types/graph_iterator_vertex.h>
#include <seqan/graph_types/graph_iterator_outedge.h>
#include <seqan/graph_types/graph_iterator_adjacency.h>
#include <seqan/graph_types/graph_iterator_edge.h>

// Graph property maps
#include <seqan/graph_types/graph_property.h>

// Specializations
#include <seqan/graph_types/graph_impl_oracle.h>	// Oracle
#include <seqan/graph_types/graph_impl_trie.h>		// Trie

// Specialized iterators
#include <seqan/graph_types/graph_iterator_bfs.h>
#include <seqan/graph_types/graph_iterator_dfs.h>

// Graph drawing and some file parsing
#include <seqan/graph_types/graph_drawing.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/graph_types/graph_utility_parsing.h>

#endif //#ifndef SEQAN_HEADER_...
