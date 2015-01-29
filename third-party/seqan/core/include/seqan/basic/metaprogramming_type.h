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
// Metaprogramming for querying and modifying types.
//
// This header contains metafunctions for querying information about types
// and modifying it, such as querying for const-ness.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_H_

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
// Metafunction IsSameType
// ----------------------------------------------------------------------------

/**
.Metafunction.IsSameType
..cat:Metaprogramming
..summary:Metaprogramming type comparison.
..signature:IsSameType<T1, T2>::Type
..signature:IsSameType<T1, T2>::VALUE
..param.T1:Left-hand argument.
..param.T2:Right-hand argument.
..returns:@Tag.Logical Values.tag.True@/$true$ if $T1$ is the same as $T2$, otherwise @Tag.Logical Values.tag.False@/$false$.
..include:seqan/basic.h
*/

template <typename Type1, typename Type2>
struct IsSameType : False {};

template <typename Type1>
struct IsSameType<Type1, Type1> : True {};

// ----------------------------------------------------------------------------
// Metafunction MakeUnsigned
// ----------------------------------------------------------------------------

/**
.Metafunction.MakeUnsigned:
..cat:Basic
..summary:Converts an integral value into an unsigned integral value.
..signature:MakeUnsigned<T>::Type
..param.T:Input integral type.
..returns.param.Type:A type without a sign of the same domain, e.g. $unsigned int$ for $T$ = $int$.
...default:$T$
..include:seqan/basic.h
 */

template <typename T>
struct MakeUnsigned
{
	typedef
		typename If<typename IsSameType<T, __int8>::Type,       __uint8,
		typename If<typename IsSameType<T, char>::Type,         unsigned char,
		typename If<typename IsSameType<T, signed char>::Type,  unsigned char,
		typename If<typename IsSameType<T, signed short>::Type, unsigned short,
		typename If<typename IsSameType<T, signed int>::Type,   unsigned int,
		typename If<typename IsSameType<T, signed long>::Type,  unsigned long,
		typename If<typename IsSameType<T, __int64>::Type,      __uint64, T
		>::Type>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeUnsigned<T const>
{
	typedef typename MakeUnsigned<T>::Type const Type;
};

// TODO(holtgrew): Internal metafunction unnecessary now?
/**
.Internal.MakeUnsigned_:
..signature:MakeUnsigned_<T>
..status:deprecated, please use @Metafunction.MakeUnsigned@
..returns:$unsigned t$ if $T$ is not $unsigned t$, otherwise $T$.
*/
template <typename T>
struct MakeUnsigned_ : MakeUnsigned<T> {};

// ----------------------------------------------------------------------------
// Metafunction MakeSigned
// ----------------------------------------------------------------------------

/**
.Metafunction.MakeSigned:
..cat:Basic
..summary:Converts an integral value into a signed integral value.
..signature:MakeSigned<T>::Type
..param.T:Input integral type.
..returns.param.Type:A type with a sign of the same domain, e.g. $int$ for $T$ = $unsigned int$.
...default:$T$
..include:seqan/basic.h
..see:Metafunction.MakeUnsigned
 */

template <typename T>
struct MakeSigned
{
	typedef
		typename If<typename IsSameType<T, char>::Type,           signed char,
		typename If<typename IsSameType<T, __int8>::Type,         __int8,
		typename If<typename IsSameType<T, unsigned char>::Type,  signed char,
		typename If<typename IsSameType<T, unsigned short>::Type, signed short,
		typename If<typename IsSameType<T, unsigned int>::Type,   signed int,
		typename If<typename IsSameType<T, unsigned long>::Type,  signed long,
		typename If<typename IsSameType<T, __uint64>::Type,       __int64, T
		>::Type>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeSigned<T const>
{
	typedef typename MakeSigned<T>::Type const Type;
};

// TODO(holtgrew): Internal metafunction unnecessary now?
/**
.Internal.MakeSigned_:
..signature:MakeSigned_<T>
..status:deprecated, please use @Metafunction.MakeSigned@
..returns:$signed t$ if $T$ is not $signed t$, otherwise $T$.
*/
template <typename T>
struct MakeSigned_ : MakeSigned<T> {};

// ----------------------------------------------------------------------------
// Metafunction RemoveReference
// ----------------------------------------------------------------------------

/**
.Metafunction.RemoveReference:
..cat:Basic
..summary:Converts a (reference) type into the same type without reference.
..signature:RemoveReference<T>::Type
..param.T:Input type.
..returns.param.Type:A corresponding non-reference type, e.g. $int$ for $T$ = $&int$.
...default:$T$
..include:seqan/basic.h
..see:Metafunction.RemoveConst
*/

// TODO(holtgrew): Internal metafunction superflous?
/**
.Internal.RemoveReference_:
..signature:RemoveReference_<T>
..status:deprecated, please use @Metafunction.RemoveReference@
..returns:$t$ if $T$ is $t &$, otherwise $T$.
*/

template <typename T>
struct RemoveReference
{
	typedef T Type;
};

template <typename T>
struct RemoveReference<T &> : RemoveReference<T> {};

template <typename T>
struct RemoveReference_ : RemoveReference<T> {};

// ----------------------------------------------------------------------------
// Metafunction RemoveConst
// ----------------------------------------------------------------------------

/**
.Metafunction.RemoveConst:
..cat:Basic
..summary:Converts a (const) type into the corresponding non-const type.
..signature:RemoveConst<T>::Type
..param.T:Input type.
..returns.param.Type:A corresponding non-const type, e.g. $int$ for $T$ = $const int$.
...default:$T$
..include:seqan/basic.h
*/

/**
.Internal.RemoveConst_:
..signature:RemoveConst_<T>
..status:deprecated, please use @Metafunction.RemoveConst@
..returns:$t$ if $T$ is $t const$, otherwise $T$.
*/

template <typename T>
struct RemoveConst
{
	typedef T Type;
};

template <typename T>
struct RemoveConst<T const> : public RemoveConst<T> {};

template <typename T>
struct RemoveConst<T &>
{
	typedef typename RemoveConst<T>::Type & Type;
};

// TODO(holtgrew): We also need a "remove inner const" meta function.
/*
template <typename T>
struct RemoveConst<T *>
{
	typedef typename RemoveConst<T>::Type * Type;
};

template <typename T, size_t I>
struct RemoveConst<T const [I]>
{
	typedef T * Type;
};
*/

// TODO(holtgrew): Internal metafunction superflous?
template <typename T>
struct RemoveConst_ : RemoveConst<T> {};

// ----------------------------------------------------------------------------
// Metafunction CopyConst_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, document.

// returns TTo const, if TFrom is const, TTo otherwise

template <typename TFrom, typename TTo>
struct CopyConst_
{
	typedef TTo Type;
};

template <typename TFrom, typename TTo>
struct CopyConst_<TFrom const, TTo>
{
	typedef TTo const Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsConst_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation.

/**
.Internal.IsConst_:
..signature:IsConst_<T>
..returns:@Tag.Logical Values.tag.True@ if $T$ is $t const$, otherwise @Tag.Logical Values.tag.False@.
*/

template <typename T>
struct IsConst_ : False
{};

template <typename T>
struct IsConst_<T const> : True
{};

// ----------------------------------------------------------------------------
// Metafunction ClassIdentifier_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation or deletion candidate.

/**
.Internal.ClassIdentifier_:
..signature:void * ClassIdentifier_<T>::getID()
..returns:A void * that identifies $T$.
...text:The returned values of two calls of $getID$ are equal if and only if
the used type $T$ was the same.
 */

template <typename T>
struct ClassIdentifier_
{
	static inline void *
	getID()
	{
		static bool _id_dummy;
		return &_id_dummy;
	}
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_METAPROGRAMMING_TYPE_H_
