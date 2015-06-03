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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// generic SIMD interface for SSE4/AVX
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_

#ifdef __SSE4_1__
//#include <immintrin.h>
#else
// SSE4.1 or greater required
// #warning "SSE4.1 instruction set not enabled"
#endif


namespace seqan {

// ============================================================================
// Useful Macros
// ============================================================================

#define SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector)                                                 \
template <typename TPosition>                                                                           \
inline typename Value<TSimdVector>::Type                                                                \
getValue(TSimdVector &vector, TPosition pos)                                                            \
{                                                                                                       \
/*                                                                                                      \
    typedef typename Value<TSimdVector>::Type TValue;                                                   \
    TValue val = (reinterpret_cast<TValue*>(&vector))[pos];                                    \
    return val;                                                                                         \
*/                                                                                                      \
    return vector[pos];                                                                                 \
}


#define SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                    \
template <typename TPosition>                                                                           \
inline typename Value<TSimdVector>::Type                                                                \
value(TSimdVector &vector, TPosition pos)                                                               \
{                                                                                                       \
    return getValue(vector, pos);                                                                       \
}

#define SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                              \
template <typename TPosition, typename TValue2>                                                         \
inline void                                                                                             \
assignValue(TSimdVector &vector, TPosition pos, TValue2 value)                                          \
{                                                                                                       \
/*                                                                                                      \
    typedef typename Value<TSimdVector>::Type TValue;                                                   \
    (reinterpret_cast<TValue*>(&vector))[pos] = value;                                                  \
*/                                                                                                      \
    vector[pos] = value;                                                                                \
}


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// define a concept and its models
// they allow us to define generic vector functions
SEQAN_CONCEPT(SimdVectorConcept, (T)) {};

// a metafunction returning the biggest supported SIMD vector
template <typename TValue>
struct SimdVector;

#define SEQAN_DEFINE_SIMD_VECTOR_(TSimdVector, TValue, TByteSize)                                       \
    typedef TValue TSimdVector __attribute__ ((__vector_size__ (TByteSize)));                           \
    template <> struct SimdVector<TValue>           { typedef TSimdVector Type; };                      \
    template <> struct Value<TSimdVector>           { typedef TValue Type; };                           \
    template <> struct LENGTH<TSimdVector>          { enum { VALUE = TByteSize / sizeof(TValue) }; };   \
    template <> struct Value<TSimdVector const>:  public Value<TSimdVector> {};                         \
    template <> struct LENGTH<TSimdVector const>: public LENGTH<TSimdVector> {};                        \
    SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector const)                                               \
    SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                        \
    SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector const)                                                  \
    SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                                  \
    template <> SEQAN_CONCEPT_IMPL((SimdVectorConcept), TSimdVector);                                   \
    template <> SEQAN_CONCEPT_IMPL((SimdVectorConcept), TSimdVector const)

#ifdef __SSE4_1__

#ifdef __AVX__
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Char,     char,           32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32SChar,    signed char,    32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32UChar,    unsigned char,  32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Short,    short,          32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UShort,   unsigned short, 32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Int,       int,            32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UInt,      unsigned int,   32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int64,     __int64,        32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt64,    __uint64,       32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Float,     float,          32);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Double,    double,         32);
#else
#ifdef __SSE3__
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Char,     char,           16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16SChar,    signed char,    16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UChar,    unsigned char,  16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Short,     short,          16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UShort,    unsigned short, 16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int,       int,            16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt,      unsigned int,   16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Int64,     __int64,        16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2UInt64,    __uint64,       16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Float,     float,          16);
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Double,    double,         16);
#endif
#endif


// ============================================================================
// Functions
// ============================================================================

#ifdef __AVX__
inline SimdVector32Char&    fill(SimdVector32Char &vector,   char x)            { return vector = _mm256_set1_epi8(x); }
inline SimdVector32SChar&   fill(SimdVector32SChar &vector,  signed char x)     { return vector = _mm256_set1_epi8(x); }
inline SimdVector32UChar&   fill(SimdVector32UChar &vector,  unsigned char x)   { return vector = _mm256_set1_epi8(x); }
inline SimdVector16Short&   fill(SimdVector16Short &vector,  short x)           { return vector = _mm256_set1_epi16(x); }
inline SimdVector16UShort&  fill(SimdVector16UShort &vector, unsigned short x)  { return vector = _mm256_set1_epi16(x); }
inline SimdVector8Int&      fill(SimdVector8Int &vector,     int x)             { return vector = _mm256_set1_epi32(x); }
inline SimdVector8UInt&     fill(SimdVector8UInt &vector,    unsigned int x)    { return vector = _mm256_set1_epi32(x); }
inline SimdVector4Int64&    fill(SimdVector4Int64 &vector,   __int64 x)         { return vector = _mm256_set1_epi64x(x); }
inline SimdVector4UInt64&   fill(SimdVector4UInt64 &vector,  __uint64 x)        { return vector = _mm256_set1_epi64x(x); }
inline SimdVector8Float&    fill(SimdVector8Float &vector,   float x)           { return vector = _mm256_set1_ps(x); }
inline SimdVector4Double&   fill(SimdVector4Double &vector,  double x)          { return vector = _mm256_set1_pd(x); }

inline void clear(SimdVector32Char &vector)     { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector32SChar &vector)    { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector32UChar &vector)    { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector16Short &vector)    { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector16UShort &vector)   { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector8Int &vector)       { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector8UInt &vector)      { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector4Int64 &vector)     { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector4UInt64 &vector)    { vector = _mm256_setzero_si256(); }
inline void clear(SimdVector8Float &vector)     { vector = _mm256_setzero_ps(); }
inline void clear(SimdVector4Double &vector)    { vector = _mm256_setzero_pd(); }

#ifdef __AVX2__
inline SimdVector32Char  shuffleVector(SimdVector32Char  const &vector, SimdVector32Char  const &indices) { return _mm256_shuffle_epi8(vector, indices); }
inline SimdVector32SChar shuffleVector(SimdVector32SChar const &vector, SimdVector32SChar const &indices) { return _mm256_shuffle_epi8(vector, indices); }
inline SimdVector32UChar shuffleVector(SimdVector32UChar const &vector, SimdVector32UChar const &indices) { return _mm256_shuffle_epi8(vector, indices); }

inline SimdVector32Char   shiftRightLogical(SimdVector32Char   const &vector, const int imm) { return _mm256_srli_epi16(vector, imm) & _mm256_set1_epi8(0xff >> imm); }
inline SimdVector32SChar  shiftRightLogical(SimdVector32SChar  const &vector, const int imm) { return _mm256_srli_epi16(vector, imm) & _mm256_set1_epi8(0xff >> imm); }
inline SimdVector32UChar  shiftRightLogical(SimdVector32UChar  const &vector, const int imm) { return _mm256_srli_epi16(vector, imm) & _mm256_set1_epi8(0xff >> imm); }
inline SimdVector16Short  shiftRightLogical(SimdVector16Short  const &vector, const int imm) { return _mm256_srli_epi16(vector, imm); }
inline SimdVector16UShort shiftRightLogical(SimdVector16UShort const &vector, const int imm) { return _mm256_srli_epi16(vector, imm); }
inline SimdVector8Int     shiftRightLogical(SimdVector8Int     const &vector, const int imm) { return _mm256_srli_epi32(vector, imm); }
inline SimdVector8UInt    shiftRightLogical(SimdVector8UInt    const &vector, const int imm) { return _mm256_srli_epi32(vector, imm); }
inline SimdVector4Int64   shiftRightLogical(SimdVector4Int64   const &vector, const int imm) { return _mm256_srli_epi64(vector, imm); }
inline SimdVector4UInt64  shiftRightLogical(SimdVector4UInt64  const &vector, const int imm) { return _mm256_srli_epi64(vector, imm); }
#else
inline SimdVector32Char  shuffleVector(SimdVector32Char  const &vector, SimdVector32Char  const &indices)
{
    return _mm256_permute2f128_si256(
        _mm256_castsi128_si256 (_mm_shuffle_epi8(_mm256_castsi256_si128(vector), _mm256_castsi256_si128(indices))),
        _mm256_castsi128_si256 (_mm_shuffle_epi8(_mm256_castsi256_si128(vector), _mm256_extractf128_si256(indices, 1))),
        0x20);
}

inline SimdVector32Char   shiftRightLogical(SimdVector32Char   const &vector, const int imm)
{
    return _mm256_permute2f128_si256(
        _mm256_castsi128_si256 (_mm_srli_epi16(_mm256_castsi256_si128(vector), imm)),
        _mm256_castsi128_si256 (_mm_srli_epi16(_mm256_extractf128_si256(vector, 1), imm)),
        0x20) & _mm256_set1_epi8(0xff >> imm);
}

#endif


template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllZeros(TSimdVector const &vector, TSimdVector const &mask)
{
#ifdef __AVX2__
    return _mm256_testz_si256(vector, mask);
#else
    return
        _mm_testz_si128(_mm256_castsi256_si128(vector), _mm256_castsi256_si128(mask)) &
        _mm_testz_si128(_mm256_extractf128_si256(vector, 1), _mm256_extractf128_si256(mask, 1));
#endif
}

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllOnes(TSimdVector const &vector)
{
#ifdef __AVX2__
    return _mm256_testc_si256(vector, _mm256_cmpeq_epi32(vector, vector));
#else
    return
        _mm_test_all_ones(_mm256_castsi256_si128(vector)) &
        _mm_test_all_ones(_mm256_extractf128_si256(vector, 1));
#endif
}

#else
#ifdef __SSE3__
inline void fill(SimdVector16Char &vector,  char x)             { vector = _mm_set1_epi8(x); }
inline void fill(SimdVector16SChar &vector, signed char x)      { vector = _mm_set1_epi8(x); }
inline void fill(SimdVector16UChar &vector, unsigned char x)    { vector = _mm_set1_epi8(x); }
inline void fill(SimdVector8Short &vector,  short x)            { vector = _mm_set1_epi16(x); }
inline void fill(SimdVector8UShort &vector, unsigned short x)   { vector = _mm_set1_epi16(x); }
inline void fill(SimdVector4Int &vector,    int x)              { vector = _mm_set1_epi32(x); }
inline void fill(SimdVector4UInt &vector,   unsigned int x)     { vector = _mm_set1_epi32(x); }
inline void fill(SimdVector2Int64 &vector,  __int64 x)          { vector = _mm_set1_epi64x(x); }
inline void fill(SimdVector2UInt64 &vector, __uint64 x)         { vector = _mm_set1_epi64x(x); }
inline void fill(SimdVector4Float &vector,   float x)           { vector = _mm_set1_ps(x); }
inline void fill(SimdVector2Double &vector,  double x)          { vector = _mm_set1_pd(x); }

inline void clear(SimdVector16Char &vector)     { vector = _mm_setzero_si128(); }
inline void clear(SimdVector16SChar &vector)    { vector = _mm_setzero_si128(); }
inline void clear(SimdVector16UChar &vector)    { vector = _mm_setzero_si128(); }
inline void clear(SimdVector8Short &vector)     { vector = _mm_setzero_si128(); }
inline void clear(SimdVector8UShort &vector)    { vector = _mm_setzero_si128(); }
inline void clear(SimdVector4Int &vector)       { vector = _mm_setzero_si128(); }
inline void clear(SimdVector4UInt &vector)      { vector = _mm_setzero_si128(); }
inline void clear(SimdVector2Int64 &vector)     { vector = _mm_setzero_si128(); }
inline void clear(SimdVector2UInt64 &vector)    { vector = _mm_setzero_si128(); }
inline void clear(SimdVector4Float &vector)     { vector = _mm_setzero_ps(); }
inline void clear(SimdVector2Double &vector)    { vector = _mm_setzero_pd(); }

inline SimdVector16Char  shuffleVector(SimdVector16Char  const &vector, SimdVector16Char  const &indices) { return _mm_shuffle_epi8(vector, indices); }
inline SimdVector16SChar shuffleVector(SimdVector16SChar const &vector, SimdVector16SChar const &indices) { return _mm_shuffle_epi8(vector, indices); }
inline SimdVector16UChar shuffleVector(SimdVector16UChar const &vector, SimdVector16UChar const &indices) { return _mm_shuffle_epi8(vector, indices); }

inline SimdVector16Char  shiftRightLogical(SimdVector16Char  const &vector, const int imm) { return _mm_srli_epi16(vector, imm) & _mm_set1_epi8(0xff >> imm); }
inline SimdVector16SChar shiftRightLogical(SimdVector16SChar const &vector, const int imm) { return _mm_srli_epi16(vector, imm) & _mm_set1_epi8(0xff >> imm); }
inline SimdVector16UChar shiftRightLogical(SimdVector16UChar const &vector, const int imm) { return _mm_srli_epi16(vector, imm) & _mm_set1_epi8(0xff >> imm); }
inline SimdVector8Short  shiftRightLogical(SimdVector8Short  const &vector, const int imm) { return _mm_srli_epi16(vector, imm); }
inline SimdVector8UShort shiftRightLogical(SimdVector8UShort const &vector, const int imm) { return _mm_srli_epi16(vector, imm); }
inline SimdVector4Int    shiftRightLogical(SimdVector4Int    const &vector, const int imm) { return _mm_srli_epi32(vector, imm); }
inline SimdVector4UInt   shiftRightLogical(SimdVector4UInt   const &vector, const int imm) { return _mm_srli_epi32(vector, imm); }
inline SimdVector2Int64  shiftRightLogical(SimdVector2Int64  const &vector, const int imm) { return _mm_srli_epi64(vector, imm); }
inline SimdVector2UInt64 shiftRightLogical(SimdVector2UInt64 const &vector, const int imm) { return _mm_srli_epi64(vector, imm); }

#ifdef __SSE4_1__
template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllZeros(TSimdVector const &vector, TSimdVector const &mask)
{
    return _mm_testz_si128(vector, mask);
}

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllOnes(TSimdVector const &vector)
{
    return _mm_test_all_ones(vector);
}

#endif
#endif
#endif


#ifdef __SSE4_1__
template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllZeros(TSimdVector const &vector)
{
    return testAllZeros(vector, vector);
}
#endif
#endif

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    std::ostream &)
inline print(std::ostream &stream, TSimdVector const &vector)
{
    stream << '<';
    for (int i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
        stream << '\t' << (unsigned)vector[i];
    stream << "\t>";
    return stream;
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_
