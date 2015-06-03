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
// This header forward declares / contains prototypes of fundamental
// metafunctions.  The rule of thumb for a metafunction to get promoted to
// "fundamental" is if it corresponds to a C++98/C++11 global typedef or is
// defined for many types.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_

#include <seqan/basic/basic_metaprogramming.h>

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
// Metafunction Value
// ----------------------------------------------------------------------------

/*!
 * @mfn Value
 * @headerfile <seqan/basic.h>
 * @brief Type of the items in the container or behind an iterator.
 *
 * @signature Value<T[, I]>::Type;
 *
 * @tparam T The type to query for its value type.
 * @tparam I Optional int, for types with multiple entries.  Defaults to 0.
 *
 * The value type of a container T is the type of the elements in T.  for example, the value type of a sequence of int
 * is int.
 */

template <typename T, const int I = 0>
struct Value;

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

/*!
 * @mfn GetValue
 * @headerfile <seqan/basic.h>
 * @brief Type for reading values.
 *
 * @signature GetValue<T>::Type;
 *
 * @tparam T Type of the value-holding object.
 *
 * Depending on T, the GetValue-type can either be Value&lt;T&gt;::Type &amp; or Value&lt;T&gt;::Type.
 *
 * @section Remarks
 *
 * <tt>GetValue</tt> is the return type of @Function.getValue@ that allows a (read-only) access to objects.  Do not
 * confuse it with value that returns a reference to the value.
 */

template <typename T>
struct GetValue;

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

/*!
 * @mfn Reference
 * @headerfile <seqan/basic.h>
 * @brief Reference type.
 *
 * @signature Reference<T>::Type;
 *
 * @tparam T A type.
 *
 * @return Type Either <tt>Value&lt;T&gt;Type &amp;</tt> or a proxy object Proxy for <tt>T</tt>.
 */

template <typename T>
struct Reference;

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn Size
 * @headerfile <seqan/basic.h>
 * @brief Type of an object that is suitable to hold size information.
 *
 * @signature Size<T>::Type;
 *
 * @tparam Type for which the size type is determined.
 *
 * @returns Type Size type of <tt>T</tt>.
 */

template <typename T>
struct Size;

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

/*!
 * @mfn Difference
 * @headerfile <seqan/basic.h>
 * @brief Difference type.
 *
 * @signature Difference<T>::Type;
 *
 * @tparam Type for which the difference type is determined.
 *
 * @returns Type Difference type of <tt>T</tt>.
 */

template <typename T>
struct Difference;

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/*!
 * @mfn Position
 * @headerfile <seqan/basic.h>
 * @brief Position type.
 *
 * @signature Position<T>::Type;
 *
 * @tparam Type for which the position type is determined.
 *
 * @returns Type position type of <tt>T</tt>.
 */

template <typename T>
struct Position;

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------


/*!
 * @mfn Spec
 * @headerfile <seqan/basic.h>
 * @brief The spec of a class.
 *
 * @signature Spec<T>::Type;
 *
 * @tparam T Type for which the spec type is determined.
 *
 * @returns Type Spec type of <tt>T</tt>.
 */

// Default case for types without Spec<>::Type specialization.

template <typename T>
struct Spec
{
    typedef void Type;
};

// Case for one template argument.

template <template <typename> class T, typename TSpec>
struct Spec<T<TSpec> >
{
    typedef TSpec Type;
};

template <typename T>
struct Spec<T const> : Spec<T>
{};

// ----------------------------------------------------------------------------
// Metafunction DeepestSpec
// ----------------------------------------------------------------------------

/*!
 * @mfn DeepestSpec
 * @headerfile <seqan/basic.h>
 * @brief The deepest spec of a class with nested template arguments.
 *
 * @signature DeepestSpec<T>::Type;
 *
 * @tparam T Type for which the deepest spec type is determined.
 *
 * @returns Type Size type of <tt>T</tt>.
 */

// Default case if not specialized for T.

template <typename T>
struct DeepestSpec
{
    typedef T Type;
};

// Default case mapping "T const" to T.

template <typename T>
struct DeepestSpec<T const> : DeepestSpec<T>
{};

// TODO(holtgrew): This can be simplified once we have variadic templates.

// Recursion for 1 argument.

template <template <typename> class T,
          typename T1>
struct DeepestSpec<T<T1> >
{
    typedef typename
        IfC<IsSameType<T1, void>::VALUE,                                  // is T1 void?
            T<T1>,                                                        // yes, end of recursion
            typename DeepestSpec<typename Spec<T<T1> >::Type >::Type      // no,  recurse
        >::Type Type;
};

// Recursion for 2 arguments.

template <template <typename, typename> class T,
          typename T1, typename T2>
struct DeepestSpec<T<T1, T2> > : DeepestSpec<typename Spec<T<T1,T2> >::Type>
{};

// Recursion for 3 arguments.

template <template <typename, typename, typename> class T,
          typename T1, typename T2, typename T3>
struct DeepestSpec<T<T1, T2, T3> >:
    DeepestSpec<typename Spec<T<T1, T2, T3> >::Type> {};

// Recursion for 4 arguments.

template <template <typename, typename, typename, typename> class T,
          typename T1, typename T2, typename T3, typename T4>
struct DeepestSpec<T<T1, T2, T3, T4> >:
    DeepestSpec<typename Spec<T<T1, T2, T3, T4> >::Type> {};

// Recursion for 5 arguments.

template <template <typename, typename, typename, typename, typename> class T,
          typename T1, typename T2, typename T3, typename T4, typename T5>
struct DeepestSpec<T<T1, T2, T3, T4, T5> > :
    DeepestSpec<typename Spec<T<T1, T2, T3, T4, T5> >::Type>
{};

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

template <typename T>
struct LENGTH;

// ----------------------------------------------------------------------------
// Metafunction MinLength
// ----------------------------------------------------------------------------

template <typename T>
struct MinLength
{
    static const unsigned VALUE = 0;
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_
