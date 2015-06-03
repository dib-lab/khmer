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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_CHUNKING_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_BASIC_CHUNKING_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// --------------------------------------------------------------------------
// Direction Tags
// --------------------------------------------------------------------------

/*!
 * @defgroup DirectionTags Direction Tags
 * @brief Tags for selecting a direction.
 */

/*!
 * @tag DirectionTags#Input
 * @headerfile <seqan/basic.h>
 * @brief Tag for selecting the input direction.
 *
 * @signature typedef Tag<Input_> Input;
 */

struct Input_;
typedef Tag<Input_> Input;

/*!
 * @tag DirectionTags#Output
 * @headerfile <seqan/basic.h>
 * @brief Tag for selecting the output direction.
 *
 * @signature typedef Tag<Output_> Output;
 */

struct Output_;
typedef Tag<Output_> Output;

/*!
 * @tag DirectionTags#Bidirectional
 * @headerfile <seqan/basic.h>
 * @brief Tag for allowing both input and output.
 *
 * @signature typedef Tag<Bidirectional_> Bidirectional;
 */

struct Bidirectional_;
typedef Tag<Bidirectional_> Bidirectional;

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Chunk
// --------------------------------------------------------------------------

/*!
 * @mfn Chunk
 * @headerfile <seqan/basic.h>
 * @brief Return the chunk type for an object.
 *
 * @signature Chunk<TObject>::Type;
 *
 * @tparam TObject The object to query for its chunk type.
 * @return Type    The chunk type of <tt>TObject</tt>.
 *
 * The default result type (if not overloaded) is @link Nothing @endlink.
 */

// Chunking is only supported for selected objects
template <typename TObject>
struct Chunk
{
    typedef Nothing Type;   // disabled by default
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function reserveChunk()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize, typename TDirection>
inline void reserveChunk(TIterator &, TSize, TDirection)
{}

// ----------------------------------------------------------------------------
// Function getChunk()
// ----------------------------------------------------------------------------

template <typename TChunk, typename TIterator, typename TDirection>
inline Nothing getChunk(TChunk &, TIterator &, TDirection)
{
    return Nothing();
}

// ----------------------------------------------------------------------------
// Function advanceChunk()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize>
inline void advanceChunk(TIterator &iter, TSize size)
{
    iter += size;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_BASIC_CHUNKING_H_
