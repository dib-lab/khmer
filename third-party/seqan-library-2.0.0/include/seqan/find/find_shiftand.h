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

#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
#define SEQAN_HEADER_FIND_SHIFTAND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ShiftAnd
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class ShiftAndPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Exact string matching using bit parallelism.
 *
 * The Shift-And algorithm is applicable to search small patterns in texts using a small alphabet.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, ShiftAnd>;
 *
 * @tparam TNeedle The needle type. Types: String
 *
 * The types of the needle and the haystack have to match.
 */

struct ShiftAnd_;
typedef Tag<ShiftAnd_> ShiftAnd;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, ShiftAnd> {
//____________________________________________________________________________
public:
    typedef unsigned int TWord;
    enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

//    Holder<TNeedle> data_host;
    String<TWord> bitMasks;            // Look up table for each character in the alphabet (called B in "Navarro")
    String<TWord> prefSufMatch;        // Set of all the prefixes of needle that match a suffix of haystack (called D in "Navarro")
    TWord needleLength;                // e.g., needleLength=33 --> blockCount=2 (iff w=32 bits)
    TWord blockCount;                // #unsigned ints required to store needle

//____________________________________________________________________________

    Pattern() {}

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
        setHost(*this, ndl);
    }

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, ShiftAnd> & me, TNeedle2 const & needle)
{
    SEQAN_CHECKPOINT
    typedef unsigned int TWord;
    typedef typename Value<TNeedle>::Type TValue;

    me.needleLength = length(needle);
    if (me.needleLength < 1)
        me.blockCount = 1;
    else
        me.blockCount = (me.needleLength - 1) / BitsPerValue<TWord>::VALUE + 1;

    clear(me.bitMasks);
    resize(me.bitMasks, me.blockCount * ValueSize<TValue>::VALUE, 0, Exact());

    for (TWord j = 0; j < me.needleLength; ++j)
        me.bitMasks[
            me.blockCount * ordValue(convert<TValue>(getValue(needle, j)))
            + j / me.MACHINE_WORD_SIZE
        ] |= (TWord)1 << (j % me.MACHINE_WORD_SIZE);

//    setValue(me.data_host, needle);

    /*
    // Debug code
    std::cout << "Alphabet size: " << ValueSize<TValue>::VALUE << std::endl;
    std::cout << "Needle length: " << me.needleLength << std::endl;
    std::cout << "Block count: " << me.blockCount << std::endl;

    for(unsigned int i=0;i<ValueSize<TValue>::VALUE;++i) {
        if ((i<97) || (i>122)) continue;
        std::cout << static_cast<char>(i) << ": ";
        for(int j=0;j<me.blockCount;++j) {
            for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
                std::cout << ((me.bitMasks[me.blockCount*i+j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
            }
        }
        std::cout << std::endl;
    }
    */
}

template <typename TNeedle, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, ShiftAnd> & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________

template <typename TNeedle>
inline TNeedle
host(Pattern<TNeedle, ShiftAnd> const & pattern)
{
SEQAN_CHECKPOINT

    typedef typename Pattern<TNeedle, ShiftAnd>::TWord TWord;
    typedef typename Value<TNeedle>::Type TValue;

    TNeedle temp;
    resize(temp, pattern.needleLength, Exact());

    TValue v = TValue();
    for (unsigned i = 0; i < length(pattern.bitMasks); i += pattern.blockCount)
    {
        for (unsigned j = 0; j < pattern.needleLength; j++)
            if ((pattern.bitMasks[i + j / pattern.MACHINE_WORD_SIZE] & (TWord)1 << (j % pattern.MACHINE_WORD_SIZE)) != (TWord)0)
                temp[j] = v;
        ++v;
    }
    return temp;
}

template <typename TNeedle>
inline TNeedle
host(Pattern<TNeedle, ShiftAnd> & pattern)
{
SEQAN_CHECKPOINT
    return host(const_cast<Pattern<TNeedle, ShiftAnd> const &>(pattern));
}

//____________________________________________________________________________

template <typename TNeedle>
inline TNeedle
needle(Pattern<TNeedle, ShiftAnd> const & pattern)
{
SEQAN_CHECKPOINT
    return host(pattern);
}

template <typename TNeedle>
inline TNeedle
needle(Pattern<TNeedle, ShiftAnd> & pattern)
{
SEQAN_CHECKPOINT
    return host(const_cast<Pattern<TNeedle, ShiftAnd> const &>(pattern));
}

//____________________________________________________________________________

template <typename TNeedle>
inline void
_patternInit (Pattern<TNeedle, ShiftAnd> & me)
{
SEQAN_CHECKPOINT
    clear(me.prefSufMatch);
    resize(me.prefSufMatch, me.blockCount, 0, Exact());
}

//____________________________________________________________________________

template <typename TFinder, typename TNeedle>
inline bool
_findShiftAndSmallNeedle(TFinder & finder, Pattern<TNeedle, ShiftAnd> & me)
{
SEQAN_CHECKPOINT
    typedef typename Value<TNeedle>::Type TValue;
    typedef unsigned int TWord;
    TWord compare = (TWord)1 << (me.needleLength - 1);
    while (!atEnd(finder))
    {
        TWord pos = ordValue(convert<TValue>(getValue(finder)));
        me.prefSufMatch[0] = ((me.prefSufMatch[0] << 1) | (TWord)1) & me.bitMasks[pos];
        if ((me.prefSufMatch[0] & compare) != 0)
        {
            _setFinderEnd(finder);
            finder -= me.needleLength - 1;
            return true;
        }
        goNext(finder);
    }
    return false;
}

template <typename TFinder, typename TNeedle>
inline bool
_findShiftAndLargeNeedle(TFinder & finder, Pattern<TNeedle, ShiftAnd> & me)
{
SEQAN_CHECKPOINT
    typedef typename Value<TNeedle>::Type TValue;
    typedef unsigned int TWord;

    TWord compare = (TWord)1 << ((me.needleLength - 1) % BitsPerValue<TWord>::VALUE);
    while (!atEnd(finder))
    {
        TWord pos = ordValue(convert<TValue>(getValue(finder)));
        TWord carry = 1;
        for(TWord block = 0; block < me.blockCount; ++block)
        {
            bool newCarry = (me.prefSufMatch[block] & ((TWord)1<< (BitsPerValue<TWord>::VALUE - 1))) != 0;
            me.prefSufMatch[block] <<= 1;
            me.prefSufMatch[block] |= carry;
            carry = newCarry;
        }
        for(TWord block = 0; block < me.blockCount; ++block)
            me.prefSufMatch[block] &= me.bitMasks[me.blockCount * pos + block];
        if ((me.prefSufMatch[me.blockCount-1] & compare) != 0)
        {
            _setFinderEnd(finder);
            finder -= me.needleLength - 1;
            return true;
        }

        /*
        // Debug code
        std::cout << "   ";
        for(int j=0;j<me.blockCount;++j) {
            for(int bit_pos=0;bit_pos<BitsPerValue<unsigned int>::VALUE;++bit_pos) {
                std::cout << ((me.prefSufMatch[j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
            }
        }
        std::cout << std::endl;
        */
        goNext(finder);
    }
    return false;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, ShiftAnd> & me) {
    SEQAN_CHECKPOINT

    if (empty(finder)) {
        _patternInit(me);
        _setFinderLength(finder, me.needleLength);
        _finderSetNonEmpty(finder);
    } else
        finder += me.needleLength;

    // Fast algorithm for needles < machine word?
    if (me.blockCount == 1)
        return _findShiftAndSmallNeedle(finder, me);
    else
        return _findShiftAndLargeNeedle(finder, me);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
