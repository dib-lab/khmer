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

// TODO(holtgrew): What's the state of this file? I remember David saying that it brought no real performance improvements.

#ifndef SEQAN_HEADER_MISC_SSE2_H
#define SEQAN_HEADER_MISC_SSE2_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifdef SEQAN_USE_SSE2_WORDS
#ifdef __SSE2__

#include <emmintrin.h>

namespace SEQAN_NAMESPACE_MAIN
{

//  #define __int128 Sse2Int128

    #ifndef SEQAN_SSE2_INT128
    #define SEQAN_SSE2_INT128
    #endif

    // ATTENTION:
    // The Sse2Int128 struct must be 16-byte aligned. Some allocators
    // don't ensure the correct alignment, e.g. the STL allocators.
    // Either avoid STL classes holding Sse2Int128 structs (maps in find_pex.h)
    // or avoid the SEQAN_USE_SSE2_WORDS define above.

    // may become obsolete when 128-bit integers will be introduced
    struct Sse2Int128
    {
    public:
        union {
            __m128i  v;
            __uint64 v64[2];
            unsigned v32[4];
        }                       data;

//____________________________________________________________________________

    public:
        Sse2Int128();
        Sse2Int128(Sse2Int128 const &);
        Sse2Int128(__m128i const &);
        Sse2Int128(__int64, __int64);
        Sse2Int128(int, int, int, int);
        Sse2Int128(short, short, short, short, short, short, short, short);

        template <typename TValue>
        Sse2Int128(TValue const &);

        //____________________________________________________________________________

        template <typename TValue>
        Sse2Int128 & operator = (TValue const &);

        //____________________________________________________________________________

        operator __m128i () const;
        operator __int64 () const;
/*      operator bool () const;
        operator __uint64 ();
        operator int ();
        operator unsigned int ();
        operator short ();
        operator unsigned short ();
        operator char ();
        operator signed char ();
        operator unsigned char ();
*/  };

    template <> struct IsSimple_<__m128i> { typedef True Type; };
    template <> struct IsSimple_<Sse2Int128> { typedef True Type; };

//____________________________________________________________________________
// clear

inline void
clear(Sse2Int128 &me)
{
    me.data.v = _mm_setzero_si128();
}

//____________________________________________________________________________
// assign

inline void
assign(Sse2Int128 &me, Sse2Int128 const &other)
{
    me.data = other.data;
}

// 1x 128bit
inline void
assign(Sse2Int128 &me, __m128i const &other)
{
    me.data.v = other;
}
inline Sse2Int128::operator __m128i () const
{
    return data.v;
}

// 64bit => 128bit
inline void
assign(Sse2Int128 &me, __int64 other)
{
    me.data.v = _mm_set_epi32(0, 0, other >> 32, other);
}
inline void
assign(Sse2Int128 &me, __uint64 other)
{
    me.data.v = _mm_set_epi32(0, 0, other >> 32, other);
}
// 64bit <= 128bit
inline Sse2Int128::operator __int64 () const
{
    return data.v64[0];
}
/*inline Sse2Int128::operator __uint64 ()
{
    return data.v64[0];
}
*/
// 32bit => 128bit
inline void
assign(Sse2Int128 &me, int other)
{
    me.data.v = _mm_set_epi32(0, 0, 0, other);
}
inline void
assign(Sse2Int128 &me, unsigned int other)
{
    me.data.v = _mm_set_epi32(0, 0, 0, other);
}
// 32bit <= 128bit
/*inline Sse2Int128::operator int ()
{
    return data.v32[0];
}
inline Sse2Int128::operator unsigned int ()
{
    return data.v32[0];
}
*/
// 16bit => 128bit
inline void
assign(Sse2Int128 &me, short other)
{
    me.data.v = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, other);
}
inline void
assign(Sse2Int128 &me, unsigned short other)
{
    me.data.v = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, other);
}
// 16bit <= 128bit
/*inline Sse2Int128::operator short ()
{
    return _mm_extract_epi16(data.v, 0);
}
inline Sse2Int128::operator unsigned short ()
{
    return _mm_extract_epi16(data.v, 0);
}
*/
// 8bit => 128bit
inline void
assign(Sse2Int128 &me, char other)
{
    me.data.v = _mm_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, other);
}
inline void
assign(Sse2Int128 &me, signed char other)
{
    me.data.v = _mm_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, other);
}
inline void
assign(Sse2Int128 &me, unsigned char other)
{
    me.data.v = _mm_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, other);
}
// 8bit <= 128bit
/*inline Sse2Int128::operator char ()
{
    return _mm_extract_epi16(data.v, 0);
}
inline Sse2Int128::operator signed char ()
{
    return _mm_extract_epi16(data.v, 0);
}
inline Sse2Int128::operator unsigned char ()
{
    return _mm_extract_epi16(data.v, 0);
}

inline Sse2Int128::operator bool () const
{
    return (__int64)(Sse2Int128)_mm_or_si128(data.v, _mm_unpackhi_epi64(data.v, data.v));
}
*/
//____________________________________________________________________________
// constructors

inline Sse2Int128::Sse2Int128()
{
    clear(*this);
}

inline Sse2Int128::Sse2Int128(Sse2Int128 const &other)
{
    assign(*this, other);
}

inline Sse2Int128::Sse2Int128(__m128i const &other)
{
    assign(*this, other);
}

// 2x 64bit
inline Sse2Int128::Sse2Int128(__int64 q1, __int64 q0)
{
    data.v = _mm_set_epi32(q1 >> 32, q1, q0 >> 32, q0);
}

// 4x 32bit
inline Sse2Int128::Sse2Int128(int q3, int q2, int q1, int q0)
{
    data.v = _mm_set_epi32(q3, q2, q1, q0);
}

// 8x 16bit
inline Sse2Int128::Sse2Int128(
                   short q7, short q6, short q5, short q4,
                   short q3, short q2, short q1, short q0)
{
    data.v = _mm_set_epi16(q7, q6, q5, q4, q3, q2, q1, q0);
}

template <typename TValue>
inline Sse2Int128::Sse2Int128(TValue const &other)
{
    assign(*this, other);
}

//____________________________________________________________________________
// operator =

template <typename TValue>
inline Sse2Int128 &
Sse2Int128::operator = (TValue const &other)
{
    assign(*this, other);
    return *this;
}

//____________________________________________________________________________
// logical operators

inline Sse2Int128
operator & (Sse2Int128 const &a, Sse2Int128 const &b)
{
    return _mm_and_si128((__m128i)a, (__m128i)b);
}
inline Sse2Int128
operator &= (Sse2Int128 &a, Sse2Int128 const &b)
{
    a.data.v = _mm_and_si128((__m128i)a, (__m128i)b);
    return a;
}

inline Sse2Int128
operator | (Sse2Int128 const &a, Sse2Int128 const &b)
{
    return _mm_or_si128((__m128i)a, (__m128i)b);
}
inline Sse2Int128
operator |= (Sse2Int128 &a, Sse2Int128 const &b)
{
    a.data.v = _mm_or_si128((__m128i)a, (__m128i)b);
    return a;
}

inline Sse2Int128
operator ^ (Sse2Int128 const &a, Sse2Int128 const &b)
{
    return _mm_xor_si128((__m128i)a, (__m128i)b);
}
inline Sse2Int128
operator ^= (Sse2Int128 &a, Sse2Int128 const &b)
{
    a.data.v = _mm_xor_si128((__m128i)a, (__m128i)b);
    return a;
}

inline Sse2Int128
operator ~ (Sse2Int128 const &a)
{
    __m128i _a = (__m128i)a;
    return _mm_xor_si128(_a, _mm_cmpeq_epi32(_a, _a));
}

//____________________________________________________________________________
// shift operators

inline Sse2Int128
operator << (Sse2Int128 const &a, int n)
{
    return _mm_or_si128(
        // n <= 64
        _mm_or_si128(
            _mm_sll_epi64((__m128i)a, (Sse2Int128)n),
            _mm_srl_epi64(
                _mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
                (Sse2Int128)(64-n)
            )
        ),
        // n >= 64
        _mm_sll_epi64(
            _mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
            (Sse2Int128)(n-64)
        )
    );
}
inline Sse2Int128
operator <<= (Sse2Int128 &a, int n)
{
    a.data.v = _mm_or_si128(
        // n <= 64
        _mm_or_si128(
            _mm_sll_epi64((__m128i)a, (Sse2Int128)n),
            _mm_srl_epi64(
                _mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
                (Sse2Int128)(64-n)
            )
        ),
        // n >= 64
        _mm_sll_epi64(
            _mm_unpacklo_epi64(_mm_setzero_si128(), (__m128i)a),
            (Sse2Int128)(n-64)
        )
    );
    return a;
}
inline Sse2Int128
operator << (Sse2Int128 const &a, unsigned int n)
{
    return a << (int)n;
}
inline Sse2Int128
operator <<= (Sse2Int128 &a, unsigned int n)
{
    return a <<= (int)n;
}

inline Sse2Int128
operator >> (Sse2Int128 const &a, int n)
{
    return _mm_or_si128(
        // n <= 64
        _mm_or_si128(
            _mm_srl_epi64((__m128i)a, (Sse2Int128)n),
            _mm_sll_epi64(
                _mm_unpackhi_epi64(a.data.v, _mm_setzero_si128()),
                (Sse2Int128)(64-n)
            )
        ),
        // n >= 64
        _mm_srl_epi64(
            _mm_unpackhi_epi64((__m128i)a, _mm_setzero_si128()),
            (Sse2Int128)(n-64)
        )
    );
}
inline Sse2Int128
operator >>= (Sse2Int128 &a, int n)
{
    a.data.v = _mm_or_si128(
        // n <= 64
        _mm_or_si128(
            _mm_srl_epi64((__m128i)a, (Sse2Int128)n),
            _mm_sll_epi64(
                _mm_unpackhi_epi64((__m128i)a, _mm_setzero_si128()),
                (Sse2Int128)(64-n)
            )
        ),
        // n >= 64
        _mm_srl_epi64(
            _mm_unpackhi_epi64((__m128i)a, _mm_setzero_si128()),
            (Sse2Int128)(n-64)
        )
    );
    return a;
}
inline Sse2Int128
operator >> (Sse2Int128 const &a, unsigned int n)
{
    return a >> (int)n;
}
inline Sse2Int128
operator >>= (Sse2Int128 &a, unsigned int n)
{
    return a >>= (int)n;
}

//____________________________________________________________________________
// artihmetic operators

//template <typename T>
inline Sse2Int128
operator + (Sse2Int128 const &a, Sse2Int128 const &b)
{
    static const __m128i carry = Sse2Int128(0, 1, 0, 0);

    union {
        __uint64 v64[2];
        __m128i v;
    } _sum;

    _sum.v = _mm_add_epi64((__m128i)a, (__m128i)b);
    if (_sum.v64[0] >= a.data.v64[0])
        return _sum.v;
    else
        return _mm_add_epi64(_sum.v, carry);
}

//template <typename T>
inline Sse2Int128
operator - (Sse2Int128 const &a, Sse2Int128 const &b)
{
    static const __m128i carry = Sse2Int128(0, 1, 0, 0);

    union {
        __uint64 v64[2];
        __m128i v;
    } _diff;

    _diff.v = _mm_sub_epi64((__m128i)a, (__m128i)b);
    if (_diff.v64[0] <= a.data.v64[0])
        return _diff.v;
    else
        return _mm_sub_epi64(_diff.v, carry);
}

// TODO: operator *, /

//____________________________________________________________________________
// compares

inline bool
operator == (Sse2Int128 const &a, Sse2Int128 const &b)
{
    __m128i _e = _mm_cmpeq_epi32((__m128i)a, (__m128i)b);
    return ((__int64)(Sse2Int128)_mm_and_si128(_e, _mm_unpackhi_epi64(_e, _e))) == ~(__int64)0;
}

inline bool
operator != (Sse2Int128 const &a, Sse2Int128 const &b)
{
    __m128i _e = _mm_cmpeq_epi32((__m128i)a, (__m128i)b);
    return ((__int64)(Sse2Int128)_mm_and_si128(_e, _mm_unpackhi_epi64(_e, _e))) != ~(__int64)0;
}

// TODO: operator <, <=, >, >=

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifdef __SSE2__
#endif //#ifdef SEQAN_USE_SSE2_WORDS

#endif //#ifndef SEQAN_HEADER_...
