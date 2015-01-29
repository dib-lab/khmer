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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Definitions for piggybacking qualities in free bits of bytes.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_QUALITIES_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_QUALITIES_H_

// TODO(holtgrew): Should the documentation be here?

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction QualityValueSize
// ----------------------------------------------------------------------------

// TODO(holtgrew): Do we want a default specialization? Should it return 0?
template <typename TValue>
struct QualityValueSize
{
    enum { VALUE = ValueSize<TValue>::VALUE };
};

template <typename TValue>
struct QualityValueSize<TValue const> : QualityValueSize<TValue>
{};

// ----------------------------------------------------------------------------
// Metafunction HasQualities
// ----------------------------------------------------------------------------

template <typename TValue>
struct HasQualities
{
    enum { VALUE = false };
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignQualityValue()
// ----------------------------------------------------------------------------

// Documentation is in alphabet_concept.h.

// ----------------------------------------------------------------------------
// Function getQualityValue()
// ----------------------------------------------------------------------------

// Documentation is in alphabet_concept.h.

// ----------------------------------------------------------------------------
// Function convertQuality()
// ----------------------------------------------------------------------------

// TODO(holtgrew): This could use some thought, what about other scales?

/**
.Function.convertQuality
..cat:Alphabets
..signature:convertQuality(c, q)
..summary:Convert an integer quality value into its ASCII representation for FASTQ (Phred scale).
..param.c:Character to store the quality in.
...type:nolink:$char$
..param.q:Value of the quality to convert.
...remarks:The quality value is an integral value between 0 and 62 (inclusive).
...type:nolink:$int$
..see:Function.getQualityValue
..include:seqan/basic.h
 */

inline 
void convertQuality(char & c, int q) 
{
    c = '!' + char(q);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_QUALITIES_H_
