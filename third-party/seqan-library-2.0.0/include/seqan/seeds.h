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
// Module for two-dimensional seeding and chaining.
// ==========================================================================

#ifndef SEQAN_HEADER_SEEDS_H
#define SEQAN_HEADER_SEEDS_H

// ===========================================================================
// Preliminaries
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <list>
#include <new>

#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/score.h>
#include <seqan/align.h>
#include <seqan/map.h>
#include <seqan/modifier.h>

// ===========================================================================
// Seeds Module
// ===========================================================================

// Basic definitions
// #include <seqan/seeds2/seeds_base.h>

// Class Seed and specializations
#include <seqan/seeds/seeds_seed_base.h>
#include <seqan/seeds/seeds_seed_simple.h>
#include <seqan/seeds/seeds_seed_diagonal.h>
#include <seqan/seeds/seeds_seed_chained.h>

// Seed extension algorithms.
#include <seqan/seeds/seeds_extension.h>

// Algorithms for chaining and merging seeds.
#include <seqan/seeds/seeds_combination.h>

// Class SeedSet, specializations, iterators.
#include <seqan/seeds/basic_iter_indirect.h>
#include <seqan/seeds/seeds_seed_set_base.h>
#include <seqan/seeds/seeds_seed_set_unordered.h>

// Banded chain alignment.
#include <seqan/seeds/banded_chain_alignment_profile.h>
#include <seqan/seeds/banded_chain_alignment_scout.h>
#include <seqan/seeds/banded_chain_alignment_traceback.h>
#include <seqan/seeds/banded_chain_alignment_impl.h>
#include <seqan/seeds/banded_chain_alignment.h>

// Global chaining algorithms
#include <seqan/seeds/seeds_global_chaining.h>

#endif  // #ifndef SEQAN_HEADER_SEEDS_H
