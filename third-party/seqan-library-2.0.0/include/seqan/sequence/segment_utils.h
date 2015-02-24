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
// Utility Functions for Sequences.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEQUENCE_SEQ_UTILS_H_
#define INCLUDE_SEQAN_SEQUENCE_SEQ_UTILS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Functor LengthLess
// --------------------------------------------------------------------------
// Compare two elements (usually segments) by their length.

template <typename T1, typename T2 = T1>
struct LengthLess
{
    bool operator() (T1 const & a, T2 const & b)
    {
        return length(a) < length(b);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function endsWith()
// ----------------------------------------------------------------------------

/*!
 * @fn endsWith
 * @brief Check whether a string is a suffix of another string.
 * @headerfile <seqan/sequence.h>
 *
 * @signature bool endsWith(str, suffix);
 *
 * @param[in] str    The string to check (A @link SegmentableConcept @endlink).
 * @param[in] suffix The suffix to check for (A @link ContainerConcept @endlink).
 *
 * @return bool <tt>true</tt> if <tt>suffix</tt> is a suffix of <tt>str</tt>.
 */

template <typename TLhs, typename TRhs>
inline bool endsWith(TLhs const & lhs, TRhs const & rhs)
{
    typename Size<TLhs const>::Type lhsLen = length(lhs);
    typename Size<TRhs const>::Type rhsLen = length(rhs);

    if (lhsLen < rhsLen)
        return false;
    return suffix(lhs, lhsLen - rhsLen) == rhs;
}

// ----------------------------------------------------------------------------
// Function startsWith()
// ----------------------------------------------------------------------------

/*!
 * @fn startsWith
 * @brief Check whether a string is a prefix of another string.
 * @headerfile <seqan/sequence.h>
 *
 * @signature bool startsWith(str, prefix);
 *
 * @param[in] str    The string to check (A @link SegmentableConcept @endlink).
 * @param[in] prefix The prefix to check for (A @link ContainerConcept @endlink).
 *
 * @return bool <tt>true</tt> if <tt>prefix</tt> is a prefix of <tt>str</tt>.
 */

// TODO(weese:) this function is doing the same as isPrefix() one should be removed
template <typename TLhs, typename TRhs>
inline bool startsWith(TLhs const & lhs, TRhs const & rhs)
{
    typename Size<TRhs const>::Type rhsLen = length(rhs);

    if (length(lhs) < rhsLen)
        return false;
    return prefix(lhs, rhsLen) == rhs;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEQUENCE_SEQ_UTILS_H_
