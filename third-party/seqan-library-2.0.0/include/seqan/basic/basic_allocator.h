// ==========================================================================
//                             basic_allocator.h
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
// Facade header for basic/allocator submodule.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALLOCATOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALLOCATOR_H_

// --------------------------------------------------------------------------
// Dependencies
// --------------------------------------------------------------------------

#include <seqan/platform.h>
#include <seqan/basic/basic_fundamental.h>
#include <seqan/basic/basic_smart_pointer.h>

// --------------------------------------------------------------------------
// Sub Module Headers
// --------------------------------------------------------------------------

// The allocator interface definitions.
#include <seqan/basic/allocator_interface.h>

// The allocator specializations.
#include <seqan/basic/allocator_simple.h>
#include <seqan/basic/allocator_singlepool.h>
#include <seqan/basic/allocator_multipool.h>
#include <seqan/basic/allocator_chunkpool.h>

// Adaption from SeqAn allocator to STL allocator.
#include <seqan/basic/allocator_to_std.h>

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_ALLOCATOR_H_
