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
// Umbrella header for the parallel module.
// ==========================================================================

#ifndef SEQAN_PARALLEL_H_
#define SEQAN_PARALLEL_H_

// ============================================================================
// Prerequisites
// ============================================================================

#include <seqan/platform.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>

#ifdef PLATFORM_WINDOWS
#include <windows.h>
#else
#include <pthread.h>
#include <errno.h>
#endif

#include <seqan/system/system_critical_section.h>   // Suspendable Queue
#include <seqan/system/system_condition.h>          // Suspendable Queue

// ----------------------------------------------------------------------------
// STL
// ----------------------------------------------------------------------------
// Use MCSTL which is part of the GCC since version 4.3

#if defined(_OPENMP) && defined(PLATFORM_GCC) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 3
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#include <numeric>
#endif // PLATFORM_GCC

#ifdef SEQAN_CXX11_STL
#include <atomic>
#include <thread>
#endif

// ============================================================================
// Module Headers
// ============================================================================

// Misc.
#include <seqan/parallel/parallel_tags.h>
#include <seqan/parallel/parallel_macros.h>

// Atomic operations.
#include <seqan/parallel/parallel_atomic_primitives.h>
#include <seqan/parallel/parallel_atomic_misc.h>
#include <seqan/parallel/parallel_lock.h>

// Splitting.
#include <seqan/parallel/parallel_splitting.h>

// Parallel variants of basic algorithms
#include <seqan/parallel/parallel_algorithms.h>

// Thread-safe / lock-free container operations.
#include <seqan/parallel/parallel_sequence.h>
#include <seqan/parallel/parallel_queue.h>
#include <seqan/parallel/parallel_queue_suspendable.h>
#include <seqan/parallel/parallel_resource_pool.h>
#include <seqan/parallel/parallel_serializer.h>

#endif  // SEQAN_PARALLEL_H_
