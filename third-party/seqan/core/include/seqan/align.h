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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_H_

// TODO(holtgrew): Usage of gapped value in align module is not consistent, need proxies in many places, reference not cleanly implemented everywhere yet.
// TODO(holtgrew): The Author: tag at the top has to be corrected in the headers of this module.
// TODO(holtgrew): Anchor Gaps must be integrated completely.
// TODO(holtgrew): Local alignments & Fragments don't work nicely together at the moment, multiLocalAlignments() needs an equivalent in the new align module.
// TODO(holtgrew): Align<>, AlignCol<> need some love and documentation.
// TODO(holtgrew): Gaps need better documentation.
// TODO(holtgrew): refinement should become graph_align and also host Graph<Alignment<>>
// TODO(holtgrew): graph_msa should become align_msa, or so, see whitepaper.
// TODO(holtgrew): The documentation and Tutorial need heavy updates, clipping alignments howto can go away.

// ============================================================================
// Prerequisites
// ============================================================================

#include <seqan/basic.h>
#include <seqan/modifier.h>  // ModifiedAlphabet<>.
#include <seqan/graph_align.h>  // TODO(holtgrew): We should not have to depend on this.

// TODO(holtgrew): Why not use priority queue from STL?
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

// ============================================================================
// Support
// ============================================================================

#include <seqan/align/gapped_value_type.h>
#include <seqan/align/align_config.h>
#include <seqan/align/align_traceback.h>
#include <seqan/align/matrix_base.h>

// ============================================================================
// Gaps & Gaps Iterator Data Structures
// ============================================================================

#include <seqan/align/gaps_base.h>
#include <seqan/align/gaps_iterator_base.h>

#include <seqan/align/gaps_array.h>
#include <seqan/align/gaps_iterator_array.h>

#include <seqan/align/gap_anchor.h>
#include <seqan/align/gaps_anchor.h>
#include <seqan/align/gaps_iterator_anchor.h>

// ============================================================================
// Alignment Data Structures and Columns
// ============================================================================

#include <seqan/align/align_metafunctions.h>
#include <seqan/align/align_cols.h>
#include <seqan/align/align_base.h>

// ============================================================================
// Alignment Algorithm Implementations.
// ============================================================================

//################################################################################
// New module
//################################################################################

// The tags have to be available everywhere so we define them centrally.
#include <seqan/align/alignment_algorithm_tags.h>

// Defines all tags needed for the DP alignment.
#include <seqan/align/dp_profile.h>

// The DP Band
#include <seqan/align/dp_band.h>

// The DP Scout
#include <seqan/align/dp_scout.h>

// Stores the score value of a particular cell in the dp matrix.
// If affine gap costs are selected one cell stores the three values
// for all three dp matrices.
#include <seqan/align/dp_cell.h>
#include <seqan/align/dp_cell_linear.h>
#include <seqan/align/dp_cell_affine.h>

// Stores the actual trace segment that was detected during traceback.
// The trace segments can be adapted into any alignment representation
// form.
#include <seqan/align/dp_trace_segment.h>
#include <seqan/align/dp_traceback_adaptor.h>

// Ensures the backwards compatibility for the global interfaces of the
// alignment algorithms. Based on the called function this selects the
// correct parameters for the new alignment module.
#include <seqan/align/dp_setup.h>

// Implements the different recursion formula of the alignment algorithms.
#include <seqan/align/dp_formula.h>
#include <seqan/align/dp_formula_linear.h>
#include <seqan/align/dp_formula_affine.h>

// Defines meta informations which determine how to compute a column and a
// certain cell for different profiles.
#include <seqan/align/dp_meta_info.h>

// Actual matrix to store the values. Uses the matrix_base.h definitions
// as a host.
#include <seqan/align/dp_matrix.h>
#include <seqan/align/dp_matrix_sparse.h>

// The navigator that based on the selected profile and band chooses the
// correct way to navigate through the matrix.
#include <seqan/align/dp_matrix_navigator.h>
#include <seqan/align/dp_matrix_navigator_score_matrix.h>
#include <seqan/align/dp_matrix_navigator_score_matrix_sparse.h>
#include <seqan/align/dp_matrix_navigator_trace_matrix.h>

// The actual implementations of the traceback and the dynamic programming that
// is used by all different alignment algorithms.
#include <seqan/align/dp_traceback_impl.h>
#include <seqan/align/dp_algorithm_impl.h>

//################################################################################
// Old module
//################################################################################

// Also, we have an implementation of Hirschberg's algorithm to compute
// alignments.
#include <seqan/align/global_alignment_hirschberg_impl.h>

// The implementations of Myers' bitvector algorithm for alignments can only
// compute alignment scores.  The combination of Hirschberg's and Myers'
// algorithm is limited in the same way.
#include <seqan/align/global_alignment_myers_impl.h>
#include <seqan/align/global_alignment_myers_hirschberg_impl.h>

// Implementations of the local alignment algorithms with declumping.  We also
// use them for the localAlignment() calls and return the best local alignment
// only.
// TODO(rmaerker): Replace this with a new implementation based on the new alignment module.
#include <seqan/align/local_alignment_waterman_eggert_impl.h>
#include <seqan/align/local_alignment_banded_waterman_eggert_impl.h>

// We carry around this implementation of Smith-Waterman because it supports
// aligning into fragment strings and alignment graphs.  Eventually, it could
// go away if Waterman-Eggert supports them.
//#include <seqan/align/local_alignment_smith_waterman_impl.h>

// ============================================================================
// Alignment Algorithm Interfaces
// ============================================================================

// The front-end functions for global alignments.
#include <seqan/align/global_alignment_unbanded.h>
#include <seqan/align/global_alignment_banded.h>

// The front-end functions for local alignments.
#include <seqan/align/local_alignment_unbanded.h>
#include <seqan/align/local_alignment_banded.h>

// The front-end for enumeration of local alignments.
#include <seqan/align/local_alignment_enumeration.h>  // documentation
#include <seqan/align/local_alignment_enumeration_unbanded.h>
#include <seqan/align/local_alignment_enumeration_banded.h>

// The front-end functions for the more specialized alignment algorithms such as
// Hirschberg, Myers and Myers-Hirschberg.
#include <seqan/align/global_alignment_specialized.h>

// ============================================================================
// Operations On Alignments
// ============================================================================

#include <seqan/align/alignment_operations.h>

#endif  // SEQAN_CORE_INCLUDE_SEQAN_ALIGN_H_
