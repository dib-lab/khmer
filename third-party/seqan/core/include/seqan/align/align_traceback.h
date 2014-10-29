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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================
// AlignTraceback object for storing the alignment traceback results.
//
// The _pump* functions for converting from AlignTrace to Gaps and Fragment
// String objects are defined where Gaps / the Alignment Graph spec is
// defined.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_TRACEBACK_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_TRACEBACK_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Specialization TraceBack
// ----------------------------------------------------------------------------

// TODO(holtgrew): Mark as internal with underscore?

/*!
 * @tag TraceBack
 * @headerfile <seqan/align.h>
 * @brief Traceback value.
 *
 * @signature struct TraceBack_;
 * @signature typedef SimpleType<unsigned char, TraceBack_> TraceBack.
 *
 * @section Remarks
 *
 * The ValueSize of <tt>TraceBack</tt> is 3.  The values are defined in the following way:
 *
 * <ul>
 *   <li>0 - Diagonal Move</li>
 *   <li>1 - Horizontal Move</li>
 *   <li>2 - Vertical Move</li>
 * </ul>
 */

/**
.Spec.TraceBack:
..cat:Alphabets
..summary: Trace back values.
..general:Class.SimpleType
..signature:TraceBack
..remarks:
...text:The @Metafunction.ValueSize@ of $TraceBack$ is 3.
The values are defined in the following way: 0=Diagonal Move, 1=Horizontal Move, 2=Vertical Move
..see:Metafunction.ValueSize
..include:seqan/align.h
*/

struct TraceBack_ {};
typedef SimpleType<unsigned char, TraceBack_> TraceBack;

template <> struct ValueSize<TraceBack>
{
    typedef __uint8 Type;
    static const Type VALUE = 3;
};

template <> struct BitsPerValue<TraceBack>
{
    typedef __uint8 Type;
    static const Type VALUE = 2;
};

// ----------------------------------------------------------------------------
// Helper Class AlignTraceback
// ----------------------------------------------------------------------------

// TODO(holtgrew): Mark as internal with underscore?

/*!
 * @class AlignTraceback
 * @headerfile <seqan/align.h>
 * @brief Data structure for storing alignment traceback.
 *
 * @signature template <typename TSize>
 *            struct AlignTraceback;
 *
 * @tparam TSize Size type to use in the traceback.
 */

/*!
 * @var TSizes AlignTraceback#sizes
 * @brief The traceback lengths.
 */

/*!
 * @var TLengths AlignTraceback#tsv
 * @brief The traceback lengths.
 */

template <typename TSize>
struct AlignTraceback
{
    // Run lengths in the align matrix.
	String<TSize> sizes;
    // Trace values: 0 = diagonal, 1 = horizontal, 2 = vertical.
	String<TraceBack> tvs;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignTracePrint()
// ----------------------------------------------------------------------------

// _alignTracePrint: this function is called by various alignment algorithm to build up the alignment during traceback

template <typename TSize, typename TSequenceH, typename TSequenceV, typename TId, typename TPos, typename TTraceValue>
inline void
_alignTracePrint(AlignTraceback<TSize> & tb,
                 TSequenceH const &,
                 TSequenceV const &,
                 TId const,
                 TPos const,
                 TId const,
                 TPos const,
                 TPos const segLen,
                 TTraceValue const tv)
{
	appendValue(tb.sizes, segLen);
	appendValue(tb.tvs, tv);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_TRACEBACK_H_
