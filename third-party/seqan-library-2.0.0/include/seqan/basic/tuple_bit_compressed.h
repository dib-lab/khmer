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
// Bit-packed tuple specialization.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_TUPLE_BIT_PACKED_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_TUPLE_BIT_PACKED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TValue>
inline bool testAllZeros(TValue const & val);

template <typename TValue>
inline bool testAllOnes(TValue const & val);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BitPackedTuple
 * @extends Tuple
 * @headerfile <seqan/basic.h>
 * @brief A fixed-size tuple of values.  Memory is saved by packing bits.
 *
 * @signature template <typename TValue, unsigned SIZE>
 *            class Tuple<TValue, SIZE, BitPacked<> >;
 *
 * @tparam TValue The value type.
 * @tparam SIZE   The length of the tuple.
 *
 * The characters are stored as bit sequence in an ordinal type (char, ..., __int64).
 *
 * @section Remarks
 *
 * Only useful for small alphabets as small tuple sizes (|Sigma|^size &lt;= 2^64 as for Dna or AminoAcid m-grams).
 */

template <unsigned char SIZE>
struct BitVector_
{
    typedef typename BitVector_<SIZE + 1>::Type Type;
};

template <> struct BitVector_<8> { typedef unsigned char Type; };
template <> struct BitVector_<16> { typedef unsigned short Type; };
template <> struct BitVector_<32> { typedef unsigned int Type; };
template <> struct BitVector_<64> { typedef __uint64 Type; };
template <> struct BitVector_<255>;

// TODO(holtgrew): There is a lot of stuff defined within the class itself. A lot of it could be moved into global functions.

// bit-packed storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
template <typename TValue, unsigned SIZE>
struct Tuple<TValue, SIZE, BitPacked<> >
{
    typedef typename BitVector_<SIZE * BitsPerValue<TValue>::VALUE>::Type TBitVector;

    static const __uint64 BIT_MASK = ((1ull << (BitsPerValue<TValue>::VALUE - 1)       ) - 1ull) << 1 | 1ull;
    static const __uint64 MASK     = ((1ull << (SIZE * BitsPerValue<TValue>::VALUE - 1)) - 1ull) << 1 | 1ull;

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

    TBitVector i;

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    // TODO(holtgrew): There is the unresolved issue whether the initialize costs critical performance. Since Tuples are PODs, it should be able to initialize Strings/arrays of them with memset().
    // TODO(weese): Use static a assertion outside of the constructor here, see SEQAN_CONCEPT_ASSERT
//    SEQAN_HOST_DEVICE Tuple() : i(0)
//    {
//        SEQAN_ASSERT_LEQ(static_cast<__uint64>(BitsPerValue<TValue>::VALUE * SIZE), static_cast<__uint64>(sizeof(TBitVector) * 8));
//    }

    // -----------------------------------------------------------------------
    // Subscription Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    template <typename TPos>
    SEQAN_HOST_DEVICE inline const TValue
    operator[](TPos k) const
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return (i >> (SIZE - 1 - k) * BitsPerValue<TValue>::VALUE) & BIT_MASK;
    }

    // -----------------------------------------------------------------------
    // Assignment Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    template <unsigned size__>
    inline Tuple & operator=(Tuple<TValue, size__, BitPacked<> > const & right)
    {
        i = right.i;
        return *this;
    }

    // TODO(holtgrew): Move the following to global functions?

    template <typename TShiftSize>
    inline TBitVector operator<<=(TShiftSize shift)
    {
        return i = (i << (shift * BitsPerValue<TValue>::VALUE)) & MASK;
    }

    template <typename TShiftSize>
    inline TBitVector operator<<(TShiftSize shift) const
    {
        return (i << (shift * BitsPerValue<TValue>::VALUE)) & MASK;
    }

    template <typename TShiftSize>
    inline TBitVector operator>>=(TShiftSize shift)
    {
        return i = (i >> (shift * BitsPerValue<TValue>::VALUE));
    }

    template <typename TShiftSize>
    inline TBitVector operator>>(TShiftSize shift) const
    {
        return i >> (shift * BitsPerValue<TValue>::VALUE);
    }

    template <typename T>
    inline void operator|=(T const & t)
    {
        i |= ordValue(t);
    }

    inline TBitVector* operator&()
    {
        return &i;
    }

    inline const TBitVector* operator&() const
    {
        return &i;
    }

    // This to be inline because elements (like this tuple) of packed structs
    // can't be arguments.
    template <typename TPos, typename TValue2>
    inline TValue2
    assignValue(TPos k, TValue2 const source)
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));

        unsigned shift = ((SIZE - 1 - k) * BitsPerValue<TValue>::VALUE);
        i = (i & ~(BIT_MASK << shift)) | (TBitVector)ordValue(source) << shift;
        return source;
    }
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// -----------------------------------------------------------------------
// Function getValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TPos>
inline TValue
getValue(Tuple<TValue, SIZE, BitPacked<> > const & me,
         TPos k)
{
    SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
    SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));

    return (me.i >> (SIZE - 1 - k) * BitsPerValue<TValue>::VALUE) & me.BIT_MASK;
}

template <typename TValue, unsigned SIZE, typename TPos>
TValue
getValue(Tuple<TValue, SIZE, BitPacked<> > & me,
         TPos k)
{
    SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
    SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));

    return (me.i >> (SIZE - 1 - k) * BitsPerValue<TValue>::VALUE) & me.BIT_MASK;
}

// -----------------------------------------------------------------------
// Function assignValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TValue2, typename TPos>
inline TValue2
assignValue(Tuple<TValue, SIZE, BitPacked<> > & me,
            TPos k,
            TValue2 const source)
{
    typedef typename Tuple<TValue, SIZE, BitPacked<> >::TBitVector TBitVector;

    SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
    SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));

    unsigned shift = ((SIZE - 1 - k) * BitsPerValue<TValue>::VALUE);
    me.i = (me.i & ~(me.BIT_MASK << shift)) | (TBitVector)ordValue(source) << shift;
    return source;
}

// -----------------------------------------------------------------------
// Function setValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TValue2, typename TPos>
inline TValue2
setValue(Tuple<TValue, SIZE, BitPacked<> > & me,
         TPos k,
         TValue2 const source)
{
    return assignValue(me, k, source);
}

// -----------------------------------------------------------------------
// Function moveValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TValue2, typename TPos>
inline TValue2
moveValue(Tuple<TValue, SIZE, BitPacked<> > & me,
          TPos k,
          TValue2 const source)
{
    return assignValue(me, k, source);
}

// ----------------------------------------------------------------------------
// Function move(), set(), assign()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline void
move(Tuple<TValue, SIZE, BitPacked<> > & t1, Tuple<TValue, SIZE, BitPacked<> > & t2)
{
    t1.i = t2.i;
}

template <typename TValue, unsigned SIZE>
inline void
set(Tuple<TValue, SIZE, BitPacked<> > & t1, Tuple<TValue, SIZE, BitPacked<> > const & t2)
{
    t1.i = t2.i;
}

template <typename TValue, unsigned SIZE>
inline void
assign(Tuple<TValue, SIZE, BitPacked<> > & t1, Tuple<TValue, SIZE, BitPacked<> > const & t2)
{
    t1.i = t2.i;
}

// -----------------------------------------------------------------------
// Function shiftLeft()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline void shiftLeft(Tuple<TValue, SIZE, BitPacked<> > & me)
{
    me <<= 1;
}

// -----------------------------------------------------------------------
// Function shiftRight()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline void shiftRight(Tuple<TValue, SIZE, BitPacked<> > & me)
{
    me >>= 1;
}

// ----------------------------------------------------------------------------
// Function testAllZeros()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline bool testAllZeros(Tuple<TValue, SIZE, BitPacked<> > const & me)
{
    return testAllZeros(me.i);
}

// ----------------------------------------------------------------------------
// Function testAllOnes()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline bool testAllOnes(Tuple<TValue, SIZE, BitPacked<> > const & me)
{
    return testAllOnes(me.i);
}

// -----------------------------------------------------------------------
// Function clear()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline void clear(Tuple<TValue, SIZE, BitPacked<> > & me)
{
    me.i = 0;
}

// ----------------------------------------------------------------------------
// Function operator&()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline Tuple<TValue, SIZE, BitPacked<> >
operator&(Tuple<TValue, SIZE, BitPacked<> > const & left,
          Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    Tuple<TValue, SIZE, BitPacked<> > tmp;
    tmp.i = left.i & right.i;
    return tmp;
}

template <typename TValue, unsigned SIZE, typename T>
inline typename Tuple<TValue, SIZE, BitPacked<> >::TBitVector
operator&(Tuple<TValue, SIZE, BitPacked<> > const & left,
          T const & right)
{
    return left.i & right;
}

// ----------------------------------------------------------------------------
// Function operator|()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline Tuple<TValue, SIZE, BitPacked<> >
operator|(Tuple<TValue, SIZE, BitPacked<> > const & left,
          Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    Tuple<TValue, SIZE, BitPacked<> > tmp;
    tmp.i = left.i | right.i;
    return tmp;
}

template <typename TValue, unsigned SIZE, typename T>
inline typename Tuple<TValue, SIZE, BitPacked<> >::TBitVector
operator|(Tuple<TValue, SIZE, BitPacked<> > const & left,
          T const & right)
{
    return left.i | right;
}

// ----------------------------------------------------------------------------
// Function operator^()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline Tuple<TValue, SIZE, BitPacked<> >
operator^(Tuple<TValue, SIZE, BitPacked<> > const & left,
          Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    Tuple<TValue, SIZE, BitPacked<> > tmp;
    tmp.i = left.i ^ right.i;
    return tmp;
}

template <typename TValue, unsigned SIZE, typename T>
inline typename Tuple<TValue, SIZE, BitPacked<> >::TBitVector
operator^(Tuple<TValue, SIZE, BitPacked<> > const & left,
          T const & right)
{
    return left.i ^ right;
}

// ----------------------------------------------------------------------------
// Function operator~()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned SIZE>
inline Tuple<TValue, SIZE, BitPacked<> >
operator~(Tuple<TValue, SIZE, BitPacked<> > const & val)
{
    Tuple<TValue, SIZE, BitPacked<> > tmp;
    tmp.i = ~val.i;
    return tmp;
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline bool operator<(Tuple<TValue, SIZE, BitPacked<> > const & left,
                      Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    return left.i < right.i;
}

template <typename TValue, unsigned SIZE>
inline bool operator<(Tuple<TValue, SIZE, BitPacked<> > & left,
                      Tuple<TValue, SIZE, BitPacked<> > & right)
{
    return left.i < right.i;
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline bool operator>(Tuple<TValue, SIZE, BitPacked<> > const & left,
                      Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    return left.i > right.i;
}

template <typename TValue, unsigned SIZE>
inline bool operator>(Tuple<TValue, SIZE, BitPacked<> > & left,
                      Tuple<TValue, SIZE, BitPacked<> > & right)
{
    return left.i > right.i;
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline bool operator<=(Tuple<TValue, SIZE, BitPacked<> > const & left,
                       Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    return !operator>(left, right);
}

template <typename TValue, unsigned SIZE>
inline bool operator<=(Tuple<TValue, SIZE, BitPacked<> > & left,
                       Tuple<TValue, SIZE, BitPacked<> > & right)
{
    return !operator>(left, right);
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline bool operator>=(Tuple<TValue, SIZE, BitPacked<> > const & left,
                       Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    return !operator<(left, right);
}

template <typename TValue, unsigned SIZE>
inline bool operator>=(Tuple<TValue, SIZE, BitPacked<> > & left,
                       Tuple<TValue, SIZE, BitPacked<> > & right)
{
    return !operator<(left, right);
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline bool operator==(Tuple<TValue, SIZE, BitPacked<> > const & left,
                       Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    return left.i == right.i;
}

template <typename TValue, unsigned SIZE>
inline bool operator==(Tuple<TValue, SIZE, BitPacked<> > & left,
                       Tuple<TValue, SIZE, BitPacked<> > & right)
{
    return left.i == right.i;
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------

// Optimized version for packed tuple using just one word.
template <typename TValue, unsigned SIZE>
inline bool operator!=(Tuple<TValue, SIZE, BitPacked<> > const & left,
                       Tuple<TValue, SIZE, BitPacked<> > const & right)
{
    return !operator==(left, right);
}

template <typename TValue, unsigned SIZE>
inline bool operator!=(Tuple<TValue, SIZE, BitPacked<> > & left,
                       Tuple<TValue, SIZE, BitPacked<> > & right)
{
    return !operator==(left, right);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_TUPLE_BIT_PACKED_H_
