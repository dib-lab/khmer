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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Basic comparison code.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_FUNDAMENTAL_COMPARISON_H_
#define SEQAN_INCLUDE_SEQAN_FUNDAMENTAL_COMPARISON_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// Forwards for Metafunctions and Functions.
template <typename T> struct ValueSize;
template <typename T> SEQAN_HOST_DEVICE inline typename ValueSize<T>::Type valueSize();
// Forwards for Metafunctions and Functions.
template <typename TValue> SEQAN_HOST_DEVICE inline typename ValueSize<TValue>::Type ordValue(TValue const & c);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn CompareType
 * @headerfile <seqan/basic.h>
 * @brief Type to convert other types for comparisons.
 *
 * @signature CompareType<T1, T2>::Type;
 *
 * @tparam T2 Type of the right operand of a comparison.
 * @tparam T1 Type of the left operand of a comparison.
 *
 * @return Type The resulting type to convert other type to.
 *
 * Comparisons are for example operators like <tt>==</tt> or <tt>&lt;</tt>.
 *
 * Do not implement, implement CompareTypeImpl instead.
 *
 * Note that there is no rule that guarantees that <tt>CompareType&lt;T1, T2&gt;::Type</tt> is the same as
 * <tt>CompareType&lt;T2, T1&gt;::Type</tt>.  It is also possible, that only one of these two types is defined.
 *
 * This metafunction is used for the implementation of comparisons that involve SimpleType.
 *
 * @see CompareTypeImpl
 */

/*!
 * @mfn CompareTypeImpl
 * @headerfile <seqan/basic.h>
 * @brief Implementation of CompareType.
 *
 * @signature CompareTypeImpl<T1, T2>::Type;
 *
 * @tparam T2 Type of the right operand of a comparison.
 * @tparam T1 Type of the left operand of a comparison.
 *
 * @return Type The type to use for the comparison.
 *
 * @see CompareType
 */

// Given two types, the CompareType is a type that both types can be cast to
// and where the results are then used to compare two values.

// step 3: choose the actual comparison type
template <typename T1, typename T2>
struct CompareTypeImpl;

template <typename T>
struct CompareTypeImpl<T, T>
{
    typedef T Type;
};

// step 2: disolve all iterator proxies (see proxy_iterator.h)
template <typename T1, typename T2>
struct CompareTypeRemoveProxy:
    CompareTypeImpl<T1, T2> {};

// step 1: remove const from types
template <typename T1, typename T2>
struct CompareType:
    CompareTypeRemoveProxy<typename RemoveConst<T1>::Type, typename RemoveConst<T2>::Type> {};

// ============================================================================
// Functions
// ============================================================================

// These functions are shortcuts to provide comparisons based on the same order
// that is imposed by ordValue

template <typename TValue1, typename TValue2>
SEQAN_HOST_DEVICE inline bool ordLess(TValue1 const & left, TValue2 const & right)
{
    return ordValue(left) < ordValue(static_cast<TValue1>(right));
}

template <typename TValue1, typename TValue2>
SEQAN_HOST_DEVICE inline bool ordEqual(TValue1 const & left, TValue2 const & right)
{
    return ordValue(left) == ordValue(static_cast<TValue1>(right));
}

template <typename TValue1, typename TValue2>
SEQAN_HOST_DEVICE inline bool ordGreater(TValue1 const & left, TValue2 const & right)
{
    return ordValue(left) > ordValue(static_cast<TValue1>(right));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_FUNDAMENTAL_COMPARISON_H_
