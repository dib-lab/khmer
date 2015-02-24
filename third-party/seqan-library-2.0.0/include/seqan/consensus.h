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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_CONSENSUS_H_
#define INCLUDE_SEQAN_CONSENSUS_H_

// ==========================================================================
// Dependencies
// ==========================================================================

#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/graph_align.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/store.h>
#include <seqan/seq_io.h>

// ==========================================================================
// Old Consensus Module
// ==========================================================================

#include <seqan/consensus/consensus_base.h>
#include <seqan/consensus/consensus_score.h>
#include <seqan/consensus/consensus_realign.h>
#include <seqan/consensus/consensus_library.h>

// ==========================================================================
// New, Store-Based Module
// ==========================================================================

#include <seqan/consensus/overlap_info_computation.h>
#include <seqan/consensus/consensus_builder.h>
#include <seqan/consensus/consensus_aligner.h>
#include <seqan/consensus/consensus_aligner_interface.h>

#endif  // #ifndef INCLUDE_SEQAN_CONSENSUS_H_
