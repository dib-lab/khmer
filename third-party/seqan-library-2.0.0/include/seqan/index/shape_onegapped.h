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

#ifndef SEQAN_HEADER_SHAPE_ONEGAPPED_H
#define SEQAN_HEADER_SHAPE_ONEGAPPED_H

namespace SEQAN_NAMESPACE_MAIN
{


    //////////////////////////////////////////////////////////////////////////////
    // shape with one gap
    //////////////////////////////////////////////////////////////////////////////

/*!
 * @class OneGappedShape
 * @extends Shape
 * @headerfile <seqan/index.h>
 *
 * @brief A variable shape with one optional gap.
 *
 * @signature template <typename TValue>
 *            class Shape<TValue, OneGappedShape>;
 *
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 *
 * A OneGappedShape must be initialized first with a valid shape. To do so, call
 * @link Shape#stringToShape @endlink.
 *
 * @see GenericShape
 */

    struct OneGappedShape;

    template <typename TValue>
    class Shape<TValue, OneGappedShape>
    {
    public:
//____________________________________________________________________________

        unsigned blockLen1;
        unsigned gapLen;
        unsigned blockLen2;

        typename Value<Shape>::Type    hValue;        // current hash value

        TValue                        leftChar;
        typename Value<Shape>::Type    factor1;
        typename Value<Shape>::Type    factor2;

//____________________________________________________________________________

        Shape():
            blockLen1(0),
            gapLen(0),
            blockLen2(0),
            hValue(0),
            leftChar(0),
            factor1(0),
            factor2(0) {}

        // c'tor for ungapped shapes
        Shape(unsigned _blockLen1, unsigned _gapLen, unsigned _blockLen2):
            blockLen1(_blockLen1),
            gapLen(_gapLen),
            blockLen2(_blockLen2),
            hValue(0),
            leftChar(0)
        {
        SEQAN_CHECKPOINT
            typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;
            factor1 = _intPow((THValue)ValueSize<TValue>::VALUE, weight(*this) - 1);
            factor2 = _intPow((THValue)ValueSize<TValue>::VALUE, blockLen2);
        }

        Shape(Shape const &other):
            blockLen1(other.blockLen1),
            gapLen(other.gapLen),
            blockLen2(other.blockLen2),
            hValue(other.hValue),
            leftChar(other.leftChar),
            factor1(other.factor1),
            factor2(other.factor2) {}

        template <typename TSpec>
        Shape(Shape<TValue, TSpec> const &other) :
                hValue(0), leftChar(0)
        {
            *this = other;
        }

        template <typename TSpec>
        Shape(GappedShape<TSpec> const &other) :
                hValue(0), leftChar(0)
        {
            *this = other;
        }

        template <typename TStringValue, typename TSpec>
        Shape(String<TStringValue, TSpec> const &bitmap) :
                hValue(0), leftChar(0)
        {
            *this = bitmap;
        }

//____________________________________________________________________________


        template <typename TStringValue, typename TSpec>
        inline Shape &
        operator=(String<TStringValue, TSpec> const &bitmap)
        {
            stringToShape(*this, bitmap);
            return *this;
        }
    };

//////////////////////////////////////////////////////////////////////////////

    template <typename TValue>
    inline typename Size< Shape<TValue, OneGappedShape> >::Type
    length(Shape<TValue, OneGappedShape> const &me)
    {
    SEQAN_CHECKPOINT
        return me.blockLen1 + me.gapLen + me.blockLen2;
    }

//____________________________________________________________________________

    template <typename TValue>
    inline typename Size< Shape<TValue, OneGappedShape> >::Type
    weight(Shape<TValue, OneGappedShape> const & me)
    {
    SEQAN_CHECKPOINT
        return me.blockLen1 + me.blockLen2;
    }

//____________________________________________________________________________

    template <typename TValue, typename TIter>
    inline typename Value< Shape<TValue, OneGappedShape> >::Type
    hash(Shape<TValue, OneGappedShape> &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;
        typedef typename Size< Shape<TValue, OneGappedShape> >::Type    TSize;

        SEQAN_ASSERT_GT(me.blockLen1, 0u);

        me.hValue = ordValue(me.leftChar = *it);
        for(TSize i = 1; i < me.blockLen1; ++i) {
            goNext(it);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        goFurther(it, me.gapLen);
        for(TSize i = 0; i < me.blockLen2; ++i) {
            goNext(it);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        return me.hValue;
    }

    template <typename TValue, typename TIter>
    inline typename Value< Shape<TValue, OneGappedShape> >::Type
    hashInit(Shape<TValue, OneGappedShape> &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;
        typedef typename Size< Shape<TValue, OneGappedShape> >::Type    TSize;

        SEQAN_ASSERT_GT(me.blockLen1, 0u);

        me.leftChar = 0;
        me.hValue = ordValue(*it);
        for(TSize i = 2; i < me.blockLen1; ++i) {
            goNext(it);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        goFurther(it, me.gapLen);
        for(TSize i = 0; i < me.blockLen2; ++i) {
            goNext(it);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        return me.hValue;
    }

    template <typename TValue, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, OneGappedShape> >::Type
    hash(Shape<TValue, OneGappedShape> &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;
        TSize blockLen1 = me.blockLen1;
        TSize blockLen2 = me.blockLen2;

        SEQAN_ASSERT_GT(me.blockLen1, 0u);

        if ((TSize)length(me) > charsLeft)
        {
            if (blockLen1 > charsLeft)
            {
                blockLen1 = charsLeft;
                blockLen2 = 0;
                if (blockLen1 == 0) return me.hValue = 0;
            } else
                if (blockLen1 + (TSize)me.gapLen > charsLeft)
                    blockLen2 = 0;
                else
                    blockLen2 = charsLeft - (blockLen1 + me.gapLen);
        }

        me.hValue = ordValue(me.leftChar = *it);
        for(TSize i = 1; i < blockLen1; ++i) {
            goNext(it);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        goFurther(it, me.gapLen);
        for(TSize i = 0; i < blockLen2; ++i) {
            goNext(it);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        // fill shape with zeros
        TSize w = (TSize)weight(me);
        for(TSize i = blockLen1 + blockLen2; i < w; ++i)
            me.hValue *= ValueSize<TValue>::VALUE;
        return me.hValue;
    }

    template <typename TValue, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, OneGappedShape> >::Type
    hashUpper(Shape<TValue, OneGappedShape> &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;
        TSize blockLen1 = me.blockLen1;
        TSize blockLen2 = me.blockLen2;

        SEQAN_ASSERT_GT(me.blockLen1, 0u);

        if ((TSize)length(me) > charsLeft)
        {
            if (blockLen1 > charsLeft)
            {
                blockLen1 = charsLeft;
                blockLen2 = 0;
            } else
                if (blockLen1 + (TSize)me.gapLen > charsLeft)
                    blockLen2 = 0;
                else
                    blockLen2 = charsLeft - (blockLen1 + me.gapLen);
        }

        me.hValue = 0;
        me.leftChar = *it;
        for(TSize i = 0; i < blockLen1; ++i) {
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
            goNext(it);
        }
        goFurther(it, me.gapLen);
        for(TSize i = 0; i < blockLen2; ++i) {
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
            goNext(it);
        }
        ++me.hValue;

        // fill shape with zeros
        TSize w = (TSize)weight(me);
        for(TSize i = blockLen1 + blockLen2; i < w; ++i)
            me.hValue *= ValueSize<TValue>::VALUE;
        return me.hValue;
    }

    template <typename TValue, typename TIter>
    inline typename Value< Shape<TValue, OneGappedShape> >::Type
    hashNext(Shape<TValue, OneGappedShape> &me, TIter &_it)
    {
        //typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;
        TIter it(_it);

        SEQAN_ASSERT_GT(me.blockLen1, 0u);

        // remove leftmost character
        me.hValue -= ordValue(me.leftChar) * me.factor1;
        me.leftChar = *it;

        // shift
        me.hValue *= ValueSize<TValue>::VALUE;

        // add emerging character
        goFurther(it, me.blockLen1 - 1);
        me.hValue += ordValue((TValue)*it) * me.factor2;

        // subtract vanishing character
        goFurther(it, me.gapLen);
        me.hValue -= ordValue((TValue)*it) * me.factor2;

        // add rightmost emerging character
        goFurther(it, me.blockLen2);
        me.hValue += ordValue((TValue)*it);

        return me.hValue;
    }


//____________________________________________________________________________
    template <typename TValue, typename TShapeString>
    inline bool
    stringToShape(
        Shape<TValue, OneGappedShape> &me,
        TShapeString const &bitmap)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TShapeString const>::Type                TIter;
        typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;

        TIter it = begin(bitmap, Standard());
        TIter itEnd = end(bitmap, Standard());

        me.blockLen1 = 0;
        me.gapLen = 0;
        me.blockLen2 = 0;

        for(; it != itEnd && *it == '0' ; ++it) ;

        for(; it != itEnd && *it == '1' ; ++it)
            ++me.blockLen1;

        for(; it != itEnd && *it == '0' ; ++it)
            ++me.gapLen;

        for(; it != itEnd && *it == '1' ; ++it)
            ++me.blockLen2;

        for(; it != itEnd && *it == '0' ; ++it) ;

        me.factor1 = _intPow((THValue)ValueSize<TValue>::VALUE, weight(me) - 1);
        me.factor2 = _intPow((THValue)ValueSize<TValue>::VALUE, me.blockLen2);

        return it == itEnd && me.blockLen1 > 0;
    }

//____________________________________________________________________________

    template <typename TShapeString, typename TValue>
    inline void
    shapeToString(
        TShapeString &bitmap,
        Shape<TValue, OneGappedShape> const &me)
    {
    SEQAN_CHECKPOINT

        clear(bitmap);
        resize(bitmap, me.blockLen1, '1');
        resize(bitmap, me.blockLen1 + me.gapLen, '0');
        resize(bitmap, me.blockLen1 + me.gapLen + me.blockLen2, '1');
    }

//____________________________________________________________________________
    template <typename TValue>
    inline void
    reverse(Shape<TValue, OneGappedShape> &me)
    {
    SEQAN_CHECKPOINT
        typedef typename Value< Shape<TValue, OneGappedShape> >::Type    THValue;

        unsigned temp = me.blockLen1;
        me.blockLen1 = me.blockLen2;
        me.blockLen2 = temp;
        me.factor2 = _intPow((THValue)ValueSize<TValue>::VALUE, me.blockLen2);
    }

}    // namespace seqan

#endif
