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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Fundamental conversion code.
// ==========================================================================

// TODO(holtgrew): Move to alphabet submodule?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_CONVERSION_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_CONVERSION_H_

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

/**
.Metafunction.Convert:
..cat:Alphabets
..summary:Return type of a conversion. 
..signature:Convert<Target, Source>::Type
..param.Target:Type the object should be converted to.
..param.Source:Type of the object that should be converted to $Target$.
..returns.param.Type:Type that is returned by @Function.convert@.
...remarks:This is either $Target$ or $Target &$:
If instances of $Source$ can be re-interpreted as instances of $Target$,
than this metafunction returns a reference, otherwise it returns $Target$, 
that is @Function.convert@ returns a temporary.
..remarks:A constant instance of $Convert$ is (ab)used as tag argument of @Function.convertImpl@.
..include:seqan/basic.h
*/

template <typename TTarget, typename TSource = void>
struct Convert
{
    typedef TTarget Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function convertImpl()
// ----------------------------------------------------------------------------

/**
.Function.convertImpl
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.convert@.
..signature:Convert convertImpl(convert, source)
..param.convert:Object that specifies the conversion.
...type:Metafunction.Convert
...remarks:A constant instance of @Metafunction.Convert@ is used to specify the conversion target.
..param.source:An object that should be converted.
..returns:$source$ converted to the type specified by convert.
...metafunction:Metafunction.Convert
..remarks:This function implements @Function.convert@. 
It is recommended to use @Function.convert@ rather than $convertImpl$.
..include:seqan/basic.h
*/

// NOTE(doering): Specialize convertImpl, use convert.
// NOTE(doering): Conversion of one char into another char happes at another place.
// NOTE(doering): Conversion of sequences happens at another place.
// NOTE(doering): Can copy or reinterpret, depending on Convert::Type

template <typename TTarget, typename T, typename TSource>
inline typename Convert<TTarget, TSource>::Type
convertImpl(Convert<TTarget, T> const,
            TSource const & source)
{
    return source;
}

// ----------------------------------------------------------------------------
// Function convert()
// ----------------------------------------------------------------------------

/**
.Function.convert
..cat:Alphabets
..summary:Converts a value into another value.
..signature:Convert convert<Target>(source)
..param.Target:The type $source$ is converted to.
..param.source:An object that is converted to $Target$.
..returns:$source$ converted to $Target$.
...remarks:If $source$ can be re-interpreted as instance of $Target$, then a reference is returned.
Otherwise the function returns a temporary object. 
...metafunction:Metafunction.Convert
..remarks:This function is implemented in @Function.convertImpl@. 
Do not specialize $convert$, specialize @Function.convertImpl@ instead.
..see:Function.convertImpl
..include:seqan/basic.h
*/

template <typename TTarget, typename TSource>
inline typename Convert<TTarget, TSource>::Type
convert(TSource const & source)
{
    return convertImpl(Convert<TTarget, TSource>(), source);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_CONVERSION_H_
