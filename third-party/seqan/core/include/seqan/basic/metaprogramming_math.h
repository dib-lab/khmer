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
// ==========================================================================
// Mathematical Metafunctions.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_MATH_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_MATH_H_

#include <seqan/platform.h>

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

/*!
 * @defgroup MetaprogrammingMath Metaprogramming Math
 * @brief Metafunctions for mathematical computations.
 */

// ----------------------------------------------------------------------------
// Metafunction Log2
// ----------------------------------------------------------------------------

/*!
 * @mfn MetaprogrammingMath#Log2
 * @brief Compute ceiled logarithm to base 2 using metaprogramming.
 *
 * @signature __uint64 Log2<NUMERUS>::VALUE;
 *
 * @tparam NUMERUS <tt>__int64</tt> value to use for the numerus.
 *
 * @return __uint64 <tt>ceil(log2(NUMERUS))</tt>
 *
 * @section Example
 *
 * @snippet core/demos/basic/metaprogramming_math.cpp log2 call
 */

/**
.Metafunction.Log2
..cat:Metaprogramming
..summary:Compute ceiled logarithm to base 2 using metaprogramming.
..signature:Log2<x>::VALUE
..param.x:The value to take the logarithm of.
...type:nolink:$__int64$
..returns:$ceil(log2(x))$.
..include:seqan/basic.h
 */

template <__int64 numerus>
struct Log2
{
    static const __uint64 VALUE = Log2<(numerus + 1) / 2>::VALUE + 1; // ceil(log_2(n))
};

// Base cases.
template <> struct Log2<1> { static const __uint64 VALUE = 0; };
template <> struct Log2<0> { static const __uint64 VALUE = 0; };

// ----------------------------------------------------------------------------
// Metafunction Log2Floor
// ----------------------------------------------------------------------------

/*!
 * @mfn MetaprogrammingMath#Log2Floor
 * @brief Compute floored logarithm to base 2 using metaprogramming.
 *
 * @signature __uint64 Log2Floor<NUMERUS>::VALUE;
 *
 * @tparam NUMERUS <tt>__int64</tt> value to use for the numerus.
 *
 * @return __uint64 <tt>floor(log2(NUMERUS))</tt>
 *
 * @section Example
 *
 * @snippet core/demos/basic/metaprogramming_math.cpp log2floor call
 */

/**
.Metafunction.Log2Floor
..cat:Metaprogramming
..summary:Compute floored logarithm to base 2 using metaprogramming.
..signature:Log2<x>::VALUE
..param.x:The value to take the logarithm of.
...type:nolink:$__int64$
..returns:$floor(log2(x))$.
..include:seqan/basic.h
 */

template <__int64 numerus>
struct Log2Floor
{
    static const __uint64 VALUE = Log2Floor<numerus / 2>::VALUE + 1;  // floor(log_2(n))
};

// Base cases.
template <> struct Log2Floor<1> { static const __uint64 VALUE = 0; };
template <> struct Log2Floor<0> { static const __uint64 VALUE = 0; };

// ----------------------------------------------------------------------------
// Metafunction Power
// ----------------------------------------------------------------------------

/*!
 * @mfn MetaprogrammingMath#Power
 * @brief Compute power of a number.
 *
 * @signature __uint64 Power<BASE, EXPONENT>::VALUE;
 *
 * @tparam BASE     The base of the term (<tt>__int64</tt>).
 * @tparam EXPONENT The exponent of the term (<tt>__int64</tt>).
 *
 * @return __uint64 b<sup>e</sup
 *
 * @snippet core/demos/basic/metaprogramming_math.cpp power call
 */

/**
.Metafunction.Power
..cat:Metaprogramming
..summary:Compute power of a number.
..signature:Power<b, e>::VALUE
..param.b:The base.
...type:nolink:$__int64$
..param.e:The exponent.
...type:nolink:$__int64$
..returns:$b^e$
..include:seqan/basic.h
 */

template <__int64 base, __int64 exponent>
struct Power {
    static const __uint64 VALUE =
            Power<base, exponent / 2>::VALUE *
            Power<base, exponent - (exponent / 2)>::VALUE;
};

// Base cases.
template <__int64 base> struct Power<base, 1> { static const __uint64 VALUE = base; };
template <__int64 base> struct Power<base, 0> { static const __uint64 VALUE = 1; };

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_MATH_H_
