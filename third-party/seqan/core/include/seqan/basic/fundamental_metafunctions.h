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
// This header forward declares / contains prototypes of fundamental
// metafunctions.  The rule of thumb for a metafunction to get promoted to
// "fundamental" is if it corresponds to a C++98/C++11 global typedef or is
// defined for many types.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_

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

/**
.Metafunction.Value:
..cat:Fundamental
..summary:Type of the items in the container or behind an iterator.
..signature:Value<T[, I]>::Type
..param.T:Type for which the value type is determined.
..param.I:Index of the entry for which to retrieve the type.
...remarks:This is only used for static-sized containers and aggregates.
...type:nolink:$int$
..returns.param.Type:Value type of $T$.
..remarks.text:
The value type of a container $T$ is the type of the elements in $T$.
For example, the value type of a sequence of $int$ is $int$.
..example.code:Value<String<char> >::Type c;  // c has type char
..include:seqan/basic.h
*/

template <typename T, const int I = 0>
struct Value;

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

/**
.Metafunction.GetValue:
..cat:Fundamental
..summary:Type for reading values. 
..signature:GetValue<T>::Type
..param.T:Type of container that holds a value.
..returns.param.Type:GetValue type of $T$.
..remarks.text:Depending on $T$, the $GetValue$-type can either be $Value<T>::Type &$ or $Value<T>::Type$.
..text:$GetValue$ is the return type of @Function.getValue@ that allows a (read-only) access to objects.
Do not confuse it with @Function.value@ that returns a @Metafunction.Reference.reference@ to the value.
..see:Metafunction.Value
..see:Function.getValue
..include:seqan/basic.h
*/

template <typename T>
struct GetValue;

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

/**
.Metafunction.Reference:
..cat:Fundamental
..summary:Reference type.
..signature:Reference<T>::Type
..param.T:A Type.
..returns.param.Type:Either $Value<T>::Type &$ or a proxy object @Class.Proxy@ for $T$.
..see:Metafunction.Value
..see:Metafunction.GetValue
..include:seqan/basic.h
*/

template <typename T>
struct Reference;

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/**
.Metafunction.Size:
..cat:Fundamental
..summary:Type of an object that is suitable to hold size information.
..signature:Size<T>::Type
..param.T:Type for which the size type is determined.
..returns.param.Type:Size type of $T$.
..remarks.text:In most cases this type is $size_t$.
..include:seqan/basic.h
*/

template <typename T>
struct Size;

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

/**
.Metafunction.Difference:
..cat:Fundamental
..summary:Type of an object that stores the difference between two iterators.
..signature:Difference<T>::Type
..param.T:Type for which the difference type is determined.
...type:Class.Iter
..returns.param.Type:Difference type of $T$.
..remarks.text:In most cases this type is $ptrdiff_t$.
..see:Metafunction.Size
..include:seqan/basic.h
*/

template <typename T>
struct Difference;

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/**
.Metafunction.Position
..cat:Fundamental
..summary:Type of an object that represents a position in a container.
..signature:Position<T>::Type
..param.T:Type for which the position type is determined.
...type:Class.Iter
...type:Class.String
..returns.param.Type:Position type of $T$.
..see:Metafunction.Iterator
..include:seqan/basic.h
*/

template <typename T>
struct Position;

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

/**
.Metafunction.Spec
..cat:Fundamental
..summary:The spec of a class. 
..signature:Spec<T>::Type
..param.T:Type for which the spec is determined.
..returns.param.Type:Spec of $T$.
..remarks:
The spec of a SeqAn type is the class that is used in template subclassing to specify the specialization. 
For example, the spec of $String<char, Alloc<> >$ is $Alloc<>$.
..remarks:
There is a default specialization for this metafunction that returns $void$.
Also, there is an implementation of Metafunction Spec for templates with one argument that returns this argument.
..include:seqan/basic.h
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

/**
.Metafunction.DeepestSpec:
..cat:Fundamental
..summary:The deepest spec of a class with nested template arguments.
..signature:DeepestSpec<T>::Type
..param.T:Type for which the deepest spec is determined.
..returns.param.Type:Deepest spec of $T$.
..remarks:The spec of a SeqAn type is the innermost class that is used in nested subclassing.
For example, the deepest spec of $Iter<..., VSTree<BottomUp<Mums> > >$ is $Mums$.
..include:seqan/basic.h
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

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_
