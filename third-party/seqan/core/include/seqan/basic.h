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
// Umbrella header for the basic module.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_H_

// TODO(holtgrew): Remove the SEQAN_CHECKPOINT macro later.
#define SEQAN_CHECKPOINT

// --------------------------------------------------------------------------
// Prerequisites
// --------------------------------------------------------------------------

#include <seqan/platform.h>

// --------------------------------------------------------------------------
// Include Sub Modules
// --------------------------------------------------------------------------

// Code for debugging and testing (assertions, test system) and profiling.
#include <seqan/basic/basic_debug.h>

// C++ Metaprogramming Support Code, generally independent of SeqAn.
#include <seqan/basic/basic_metaprogramming.h>

// Fundamental meta and global functions.  This is what makes SeqAn SeqAn.
#include <seqan/basic/basic_fundamental.h>

// More advanced debug system constructs.
// TODO(holtgrew): Move into basic_debug subsystem, some stuff from metaprogramming and fundamental required, those should not depend on debug system.
#include <seqan/basic/test_system.h>

// SeqAn Concept Checking Library (ported from Boost).
#include <seqan/basic/basic_concept.h>

// Alphabet concept and biological implementations.
#include <seqan/basic/basic_alphabet.h>

// Aggregate data types (pairs, triples, tuples).
#include <seqan/basic/basic_aggregate.h>

// Memory allocation code.
#include <seqan/basic/basic_allocator.h>

// High level parallelism support.
#include <seqan/basic/basic_parallelism.h>

// Mathematical functions and utilities.
#include <seqan/basic/basic_math.h>

// Smart pointers, including Holder<> class hierarchy.
#include <seqan/basic/basic_smart_pointer.h>

// Iterator concept and implementation.
#include <seqan/basic/basic_iterator.h>

// Proxy class and supporting code.
#include <seqan/basic/basic_proxy.h>

// Container concept and supporting code.
#include <seqan/basic/basic_container.h>

// Remaining code with cyclic dependencies.
#include <seqan/basic/basic_tangle.h>

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_H_
