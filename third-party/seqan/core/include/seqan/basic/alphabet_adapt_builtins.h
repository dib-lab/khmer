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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Adaptions of builting types such as bool, int, but also "builtin-level"
// user defined types such as wchar_t, __int64, __uint64 to the alphabet
// concepts they are in.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_BASIC_ALPHABET_ADAPT_BUILTINS_H_
#define SEQAN_CORE_INCLUDE_BASIC_ALPHABET_ADAPT_BUILTINS_H_

#include <limits>

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
// Metafunctions MaxValue_, MinValue_
// ----------------------------------------------------------------------------

// We would want to have this here, however this is not possible with the
// current implementation.

// ----------------------------------------------------------------------------
// Metafunction BitsPerValue
// ----------------------------------------------------------------------------

template <>
struct BitsPerValue<bool>
{
    typedef int Type;
    enum { VALUE = 1 };
};

// ----------------------------------------------------------------------------
// Metafunction IsCharType
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should probably become a concept.

/**
.Metafunction.IsCharType
..cat:Alphabets
..summary:Return whether the argument is $char$, $wchar_t$, $char const$, or $wchar_t const$.
..signature:IsCharType<T>::Type
..signature:IsCharType<T>::VALUE
..param.T:Type to check type of.
..remarks:This metafunction is used to enable and disable templated adaptions of arrays to sequences for builtin character types only.
..remarks:The return value is $True$/$true$ for $char$, $wchar_t$, $char const$, and $wchar_t const$.
..include:seqan/sequence.h
*/

template <typename T>
struct IsCharType;

template <typename T>
struct IsCharType
{
    typedef False Type;
    enum { VALUE = 0 };
};

template <typename T>
struct IsCharType<T const>
    : IsCharType<T> {};

template <>
struct IsCharType<char>
{
    typedef True Type;
    enum { VALUE = 1 };
};

template <>
struct IsCharType<wchar_t>
{
    typedef True Type;
    enum { VALUE = 1 };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function gapValueImpl()                                               [char]
// ----------------------------------------------------------------------------

inline char const &
gapValueImpl(char *)
{
    static char const _gap = '-';
    return _gap;
}

inline char const &
gapValueImpl(char const *)
{
    static char const _gap = '-';
    return _gap;
}

// ----------------------------------------------------------------------------
// Function unknownValueImpl()                                           [char]
// ----------------------------------------------------------------------------

inline char const &
unknownValueImpl(char *)
{
    static char const _unknown = 'N';
    return _unknown;
}

inline char const &
unknownValueImpl(char const *)
{
    static char const _unknown = 'N';
    return _unknown;
}

// ----------------------------------------------------------------------------
// Function supremumValueImpl()
// ----------------------------------------------------------------------------

template <typename T>
inline T const &
supremumValueImpl(T *)
{
    static T const x = MaxValue<T>::VALUE;
    return x;
}

inline long double const &
supremumValueImpl(long double *)
{
#ifdef PLATFORM_WINDOWS
    static long double const _value = ::std::numeric_limits<long double>::infinity( );
#else
    static long double const _value = 1.7976931348623157e+308;
#endif
    return _value;
}

inline double const &
supremumValueImpl(double *)
{
#ifdef PLATFORM_WINDOWS
    static double const _value = ::std::numeric_limits<double>::infinity( );
#else
    static double const _value = 1.7976931348623157e+308;
#endif
    return _value;
}
inline float const &
supremumValueImpl(float *)
{
#ifdef PLATFORM_WINDOWS
    static float const _value = ::std::numeric_limits<float>::infinity( );
#else
    static float const _value = 3.40282347e+38F;
#endif
    return _value;
}

// ----------------------------------------------------------------------------
// Function infimumValueImpl()
// ----------------------------------------------------------------------------

template <typename T>
inline T const &
infimumValueImpl(T *)
{
    static T const x = MinValue<T>::VALUE;
    return x;
}

inline float const &
infimumValueImpl(float *)
{
#ifdef PLATFORM_WINDOWS
    static float const _value = -::std::numeric_limits<float>::infinity( );
#else
    static float const _value = -3.40282347e+38F;
#endif
    return _value;
}

inline double const &
infimumValueImpl(double *)
{
#ifdef PLATFORM_WINDOWS
    static double const _value = -::std::numeric_limits<double>::infinity( );
#else
    static double const _value = -1.7976931348623157e+308;
#endif
    return _value;
}

inline long double const &
infimumValueImpl(long double *)
{
#ifdef PLATFORM_WINDOWS
    static long double const _value = -::std::numeric_limits<long double>::infinity( );
#else
    static long double const _value = -1.7976931348623157e+308;
#endif
    return _value;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_BASIC_ALPHABET_ADAPT_BUILTINS_H_
