// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Tobias Rausch <rausch@embl.de>
// Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
// ==========================================================================
// Umbrella header for the moduel graph_align.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_GRAPH_ALIGN_H_
#define SEQAN_INCLUDE_SEQAN_GRAPH_ALIGN_H_

// External STL
#include <map>

// Seqan
#include <seqan/score.h>
#include <seqan/align/fragment.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>


// Alignment graph
#include <seqan/graph_align/graph_impl_align.h>
#include <seqan/graph_align/graph_impl_align_adapt.h>

// Interval trees
#include <seqan/misc/interval_tree.h>

// Refinement
//#include <seqan/graph_align/graph_algorithm_refine.h>
#include <seqan/graph_align/graph_algorithm_refine_scoring.h>
#include <seqan/graph_align/graph_algorithm_refine_fragment.h>
#include <seqan/graph_align/graph_algorithm_refine_aligngraph.h>
#include <seqan/graph_align/graph_algorithm_refine_align.h>
//#include <seqan/graph_align/graph_algorithm_refine_exact.h>
#include <seqan/graph_align/graph_algorithm_refine_exact_iterative.h>
#include <seqan/graph_align/graph_algorithm_refine_inexact.h>
#include <seqan/graph_align/graph_algorithm_refine_annotation.h>

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_GRAPH_ALIGN_H_
