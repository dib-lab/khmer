// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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

#ifndef SEQAN_HEADER_INDEX_H
#define SEQAN_HEADER_INDEX_H

// ==========================================================================
// Prerequisites.
// ==========================================================================

#include <seqan/sequence.h>
#include <seqan/pipe.h>
#include <seqan/modifier.h>
#include <seqan/parallel.h>

#include <seqan/find.h>
#include <seqan/misc/set.h>
#include <seqan/misc/memset.h>
#include <seqan/misc/bit_twiddling.h>

#include <climits>
#include <functional>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cmath>    // isnan, isinf
#include <string.h> // memset

#ifdef __SSE4_2__
#include <nmmintrin.h>
#endif

// ==========================================================================
// Index construction.
// ==========================================================================

#include <seqan/index/index_forwards.h>

#ifdef SEQAN_PARALLEL
#include <seqan/parallel.h>
#endif  // #ifdef SEQAN_PARALLEL

#include <seqan/index/index_base.h>

// ----------------------------------------------------------------------------
// q-Gram index creator.
// ----------------------------------------------------------------------------

#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>
#include <seqan/index/shape_onegapped.h>
#include <seqan/index/shape_predefined.h>
#include <seqan/index/shape_threshold.h>
#include <seqan/index/index_qgram.h>
#include <seqan/index/index_qgram_openaddressing.h>
//#include <seqan/index/index_qgram_nested.h>

// ----------------------------------------------------------------------------
// Suffix array creators.
// ----------------------------------------------------------------------------

#include <seqan/index/radix.h>
#include <seqan/index/index_sa_btree.h>
#include <seqan/index/index_sa_lss.h>
#include <seqan/index/index_sa_mm.h>
#include <seqan/index/index_sa_qsort.h>
#include <seqan/index/index_sa_bwtwalk.h>

#include <seqan/index/pump_extender3.h>
#include <seqan/index/pipe_merger3.h>
#include <seqan/index/index_skew3.h>

#include <seqan/index/pump_extender7.h>
#include <seqan/index/pipe_merger7.h>
#include <seqan/index/index_skew7.h>

#include <seqan/index/pump_separator7.h>
#include <seqan/index/index_skew7_multi.h>

// ----------------------------------------------------------------------------
// Enhanced table creators.
// ----------------------------------------------------------------------------

#include <seqan/index/pump_lcp_core.h>
#include <seqan/index/index_lcp.h>
#include <seqan/index/index_lcp_tree.h>

#include <seqan/index/index_childtab.h>
#include <seqan/index/index_bwt.h>

// ==========================================================================
// Index usage.
// ==========================================================================

#include <seqan/index/index_shims.h>

// ----------------------------------------------------------------------------
// Virtual string trees.
// ----------------------------------------------------------------------------

#include <seqan/index/index_esa_base.h>
#include <seqan/index/index_esa_stree.h>
#include <seqan/index/index_wotd.h>
#include <seqan/index/index_dfi.h>
#include <seqan/index/index_sa_stree.h>
#include <seqan/index/index_sa_truncated.h>
#include <seqan/index/index_qgram_stree.h>
#include <seqan/index/index_qgram_bucketrefinement.h>

// ----------------------------------------------------------------------------
// Rank dictionaries.
// ----------------------------------------------------------------------------

#include <seqan/index/index_fm_rank_dictionary_base.h>
#include <seqan/index/index_fm_rank_dictionary_naive.h>
#include <seqan/index/index_fm_rank_dictionary_levels.h>
#include <seqan/index/index_fm_right_array_binary_tree.h>
#include <seqan/index/index_fm_right_array_binary_tree_iterator.h>
#include <seqan/index/index_fm_rank_dictionary_wt.h>

// ----------------------------------------------------------------------------
// Sparse strings.
// ----------------------------------------------------------------------------

#include <seqan/index/index_fm_sparse_string.h>

// ----------------------------------------------------------------------------
// FM index.
// ----------------------------------------------------------------------------

#include <seqan/index/index_fm_lf_table.h>
#include <seqan/index/index_fm_compressed_sa.h>
#include <seqan/index/index_fm_compressed_sa_iterator.h>
#include <seqan/index/index_fm.h>
#include <seqan/index/index_fm_stree.h>

// ----------------------------------------------------------------------------
// Suffix tree algorithms.
// ----------------------------------------------------------------------------

#include <seqan/index/index_esa_algs.h>
#include <seqan/index/index_esa_algs_multi.h>
#include <seqan/index/index_esa_drawing.h>
#include <seqan/index/repeat_base.h>

// ----------------------------------------------------------------------------
// Pizza & Chili interface (compressed indices).
// ----------------------------------------------------------------------------

#include <seqan/index/index_pizzachili.h>
#include <seqan/index/index_pizzachili_find.h>

// ----------------------------------------------------------------------------
// Shawarma interface (suffix array creators).
// ----------------------------------------------------------------------------

#include <seqan/index/index_shawarma.h>

// ----------------------------------------------------------------------------
// View and device indices.
// ----------------------------------------------------------------------------

#include <seqan/index/index_view.h>
// NOTE(esiragusa): CUDA FM-index is broken.
//#ifdef PLATFORM_CUDA
//#include <seqan/index/index_device.h>
//#endif

// ==========================================================================
// Finder interface.
// ==========================================================================

// ----------------------------------------------------------------------------
// Index based finders.
// ----------------------------------------------------------------------------

#include <seqan/index/find_index_esa.h>
#include <seqan/index/find_index_approx.h>
#include <seqan/index/find_index_qgram.h>
#include <seqan/index/find_index.h>
#include <seqan/index/find_quasar.h>
#include <seqan/index/find_swift.h>
#include <seqan/index/find_pigeonhole.h>
#include <seqan/index/find_backtracking.h>

// ----------------------------------------------------------------------------
// Non-standard index based finders.
// ----------------------------------------------------------------------------

#include <seqan/index/find2_base.h>
#include <seqan/index/find2_index.h>
#include <seqan/index/find2_backtracking.h>
#include <seqan/index/find2_vstree_factory.h>
#include <seqan/index/find2_index_multi.h>
#include <seqan/index/find2_functors.h>

// ----------------------------------------------------------------------------
// Lambda interface.
// ----------------------------------------------------------------------------

#ifdef SEQAN_CXX11_STANDARD
#include <seqan/index/find_index_lambda.h>
#endif

#endif //#ifndef SEQAN_HEADER_...
