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

#ifndef SEQAN_HEADER_SHAPE_GAPPED_H
#define SEQAN_HEADER_SHAPE_GAPPED_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// HardwiredShape allows compiler-time defined gapped shape

/*!
 * @class HardwiredShape
 * @extends GappedShape
 * @headerfile <seqan/index.h>
 * @brief A structure to define a fixed gapped shape.
 *
 * @signature template <typename TValue, typename P1[, typename P2[, ..., typename Pn]]>
 *            class Shape<TValue, HardwiredShape<P1[, P2[, ..., Pn]]> >;
 *
 * @tparam TValue The type to use for the values.
 * @tparam Pi     Pi is the distance of the i'th '1' to the next '1' in the shape.At
 *                most 20 parameters are allowed, so the maximal shape weight is 21.
 *
 * You can use this structure to define your one gapped shapes in conjunction
 * with @link GappedShape @endlink.
 *
 * The shape <tt>1100101</tt> corresponds to <tt>HardwiredShape<1,3,2></tt>.
 *
 * Note The following predefined shapes are already available in <tt>seqan/index/shape_predefined.h</tt>:
 *
 * @include include/seqan/index/shape_predefined.h
 *
 * @see GappedShape
 * @see GenericShape
 */

    // Pxx = spaces between '1's
    template <
        int P00 = 0, int P01 = 0, int P02 = 0, int P03 = 0, int P04 = 0,
        int P05 = 0, int P06 = 0, int P07 = 0, int P08 = 0, int P09 = 0,
        int P10 = 0, int P11 = 0, int P12 = 0, int P13 = 0, int P14 = 0,
        int P15 = 0, int P16 = 0, int P17 = 0, int P18 = 0, int P19 = 0
    >
    struct HardwiredShape {
        static const int DIFFS[];
    };

    template <
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19
    >
    const int HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    >::DIFFS[] = {
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19, 0
    };


//////////////////////////////////////////////////////////////////////////////
// LENGTH meta-function for fixed gapped shapes

    template <>
    struct LENGTH< HardwiredShape<
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0> >
    {
        enum { VALUE = 1 };
    };

    template <
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19
    >
    struct LENGTH< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19> >
    {
        enum { VALUE = LENGTH< HardwiredShape<
            P01,P02,P03,P04,P05,
            P06,P07,P08,P09,P10,
            P11,P12,P13,P14,P15,
            P16,P17,P18,P19, 0 > >::VALUE + P00 };
    };

    template <
        typename TValue,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19
    >
    struct LENGTH< Shape<TValue, GappedShape< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    > > > >:
    LENGTH< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    > > {};


//////////////////////////////////////////////////////////////////////////////
// WEIGHT meta-function for fixed gapped shapes

    template <>
    struct WEIGHT< HardwiredShape<
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0> >
    {
        enum { VALUE = 1 };
    };

    template <
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19
    >
    struct WEIGHT< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19> >
    {
        enum { VALUE = WEIGHT< HardwiredShape<
            P01,P02,P03,P04,P05,
            P06,P07,P08,P09,P10,
            P11,P12,P13,P14,P15,
            P16,P17,P18,P19, 0 > >::VALUE + 1 };
    };

    template <
        typename TValue,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19
    >
    struct WEIGHT< Shape<TValue, GappedShape< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    > > > >:
    WEIGHT< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    > > {};


//////////////////////////////////////////////////////////////////////////////

/*!
 * @class GenericShape
 * @extends Shape
 * @headerfile <seqan/index.h>
 * @brief A variable gapped shape.
 *
 * @signature template <typename TValue>
 *            class Shape<TValue, GenericShape>;
 *
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 *
 * A GenericShape must be initialized first with a valid shape.  To do so, call
 * @link Shape#stringToShape @endlink.
 *
 * @see GappedShape
 * @see OneGappedShape
 */
    //////////////////////////////////////////////////////////////////////////////
    // variable gapped shape
    //////////////////////////////////////////////////////////////////////////////

    template <typename TValue>
    class Shape<TValue, GenericShape>
    {
    public:
//____________________________________________________________________________

        unsigned span;
        unsigned weight;
        String<int> diffs;

        typename Value<Shape>::Type    hValue;        // current hash value
//____________________________________________________________________________

/*!
 * @fn GenericShape::Shape
 *
 * @brief Constructor
 *
 * @signature Shape::Shape();
 * @signature Shape::Shape(q);
 * @signature Shape::Shape(shape);
 * @signature Shape::Shape(bitmap);
 * @signature Shape::Shape(predefined);
 *
 * @param predefined Any instance of a predefined shape spec (e.g.
 *                   <tt>ShapePatternHunter</tt>).
 * @param q Creates an ungapped q-gram.
 * @param shape Any other gapped/ungapped shape.
 * @param bitmap Bitmap string. Sequence of '0's and '1's. See @link Shape#stringToShape @endlink
 *
 * @see HardwiredShape
 */
        Shape():
            span(0),
            weight(0),
            hValue(0) {}

        // c'tor for ungapped shapes
        Shape(unsigned _span):
            span(_span),
            weight(_span),
            hValue(0)
        {
        SEQAN_CHECKPOINT
            if (_span > 0)
            {
                resize(diffs, _span - 1);
                for(unsigned i = 0; i < _span - 1; ++i)
                    diffs[i] = 1;
            }
        }

        Shape(Shape const &other):
            span(other.span),
            weight(other.weight),
            diffs(other.diffs),
            hValue(other.hValue) {}

        template <typename TSpec>
        Shape(Shape<TValue, TSpec> const &other)
        {
            *this = other;
        }

        template <typename TSpec>
        Shape(GappedShape<TSpec> const &other)
        {
            *this = other;
        }

        template <typename TStringValue, typename TSpec>
        Shape(String<TStringValue, TSpec> const &bitmap)
        {
            *this = bitmap;
        }

//____________________________________________________________________________

        template <unsigned q>
        inline Shape &
        operator=(Shape<TValue, UngappedShape<q> > const &other)
        {
            span = length(other);
            weight = seqan::weight(other);
            if (weight > 0)
            {
                resize(diffs, weight - 1);
                for(unsigned i = 0; i < weight - 1; ++i)
                    diffs[i] = 1;
            }
            else
                clear(diffs);
            hValue = other.hValue;
            return *this;
        }

        template <typename TSpec>
        inline Shape &
        operator=(Shape<TValue, GappedShape<TSpec> > const &other)
        {
            span = other.span;
            weight = other.weight;
            diffs = other.diffs;
            hValue = other.hValue;
            return *this;
        }

        template <typename TSpec>
        inline Shape &
        operator=(GappedShape<TSpec> const)
        {
            typedef Shape<TValue, GappedShape<TSpec> > TShape;
            return *this = TShape();
        }

        template <typename TStringValue, typename TSpec>
        inline Shape &
        operator=(String<TStringValue, TSpec> const &bitmap)
        {
            stringToShape(*this, bitmap);
            return *this;
        }
    };


//////////////////////////////////////////////////////////////////////////////

/*!
 * @class GappedShape
 * @extends Shape
 * @headerfile <seqan/index.h>
 *
 * @brief A fixed gapped shape.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class Shape<TValue, GappedShape<TSpec> >;
 *
 * @tparam TSpec A structure to store the shape at compile-time. Types:
 *               @link HardwiredShape @endlink
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 *
 * There are predefined shapes in <tt>index/shape_predefined.h</tt>. You can
 * simply use them with <tt>Shape<TValue, ShapePatternHunter></tt> for example.
 *
 * @see HardwiredShape
 * @see GenericShape
 */

    //////////////////////////////////////////////////////////////////////////////
    // fixed gapped shape
    //////////////////////////////////////////////////////////////////////////////

    template <typename TValue, typename TSpec>
    class Shape<TValue, GappedShape<TSpec> >
    {
    public:
//____________________________________________________________________________

        typedef GappedShape<TSpec>    TShapeSpec;

        enum { span = LENGTH<Shape>::VALUE };
        enum { weight = WEIGHT<Shape>::VALUE };
        const int *diffs;

        typename Value<Shape>::Type    hValue;        // current hash value
//____________________________________________________________________________

        Shape():
            diffs(TSpec::DIFFS),
            hValue(0) {}

        Shape(Shape const &other):
            diffs(other.diffs),
            hValue(other.hValue) {}
//____________________________________________________________________________

        inline Shape &
        operator=(Shape const &other)
        {
            hValue = other.hValue;
      return *this;
        }
    };

//////////////////////////////////////////////////////////////////////////////

    template <typename TValue, typename TSpec>
    inline typename Size< Shape<TValue, GappedShape<TSpec> > >::Type
    weight(Shape<TValue, GappedShape<TSpec> > const & me)
    {
    SEQAN_CHECKPOINT
        return me.weight;
    }

//____________________________________________________________________________

    template <typename TValue, typename TIter>
    inline typename Value< Shape<TValue, GenericShape> >::Type
    hash(Shape<TValue, GenericShape> &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, GenericShape> >::Type    THValue;
        typedef typename Size< Shape<TValue, GenericShape> >::Type    TSize;

        SEQAN_ASSERT_GT((unsigned)me.weight, 0u);

        me.hValue = ordValue((TValue)*it);
        TSize iEnd = me.weight - 1;
        for(TSize i = 0; i < iEnd; ++i) {
            goFurther(it, me.diffs[i]);
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        return me.hValue;
    }

    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, GappedShape<TSpec> > >::Type
    hash(Shape<TValue, GappedShape<TSpec> > &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, GappedShape<TSpec> > >::Type    THValue;

        if (charsLeft >= (TSize)me.span)
            return hash(me, it);

        SEQAN_ASSERT_GT((unsigned)me.weight, 0u);

        TSize i = 0;
        if (charsLeft > 0) {
            me.hValue = ordValue((TValue)*it);
            TSize d;
            while (charsLeft > (d = me.diffs[i])) {
                charsLeft -= d;
                goFurther(it, d);
                me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
                ++i;
            }
            ++i;
        } else
            return me.hValue = 0;

        // fill shape with zeros
        for(; i < (TSize)me.weight; ++i)
            me.hValue *= ValueSize<TValue>::VALUE;
        return me.hValue;
    }

    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, GappedShape<TSpec> > >::Type
    hashUpper(Shape<TValue, GappedShape<TSpec> > &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, GappedShape<TSpec> > >::Type    THValue;

        if (charsLeft >= (TSize)me.span) {
            hash(me, it);
            return ++me.hValue;
        }

        SEQAN_ASSERT_GT((unsigned)me.weight, 0u);

        TSize i = 0;
        if (charsLeft > 0) {
            me.hValue = ordValue((TValue)*it);
            TSize d;
            while (charsLeft > (d = me.diffs[i])) {
                charsLeft -= d;
                goFurther(it, d);
                me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
                ++i;
            }
            ++i;
            ++me.hValue;
        } else
            me.hValue = 1;

        // fill shape with zeros
        for(; i < (TSize)me.weight; ++i)
            me.hValue *= ValueSize<TValue>::VALUE;
        return me.hValue;
    }

//____________________________________________________________________________

    template <typename THValue, typename TValue, typename TIter>
    inline THValue
    _hashHardwiredShape(THValue hash, TIter &, TValue const, HardwiredShape<
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0 > const)
    {
        return hash;
    }

    template <
                 int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19,
        typename THValue, typename TValue, typename TIter
    >
    inline THValue
    _hashHardwiredShape(THValue hash, TIter &it, TValue const, HardwiredShape<
         1 ,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19 > const)
    {
        ++it;
        return _hashHardwiredShape(hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
            it, TValue(), HardwiredShape<
                P01,P02,P03,P04,P05,
                P06,P07,P08,P09,P10,
                P11,P12,P13,P14,P15,
                P16,P17,P18,P19, 0 >());
    }

    template <
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19,
        typename THValue, typename TValue, typename TIter
    >
    inline THValue
    _hashHardwiredShape(THValue hash, TIter &it, TValue const, HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19 > const)
    {
        it += P00;
        return _hashHardwiredShape(hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
            it, TValue(), HardwiredShape<
                P01,P02,P03,P04,P05,
                P06,P07,P08,P09,P10,
                P11,P12,P13,P14,P15,
                P16,P17,P18,P19, 0 >());
    }

    template <
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19,
        typename TValue, typename TIter
    >
    inline typename Value< Shape<TValue, GappedShape< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    > > > >::Type
    hash(Shape<TValue, GappedShape< HardwiredShape<
        P00,P01,P02,P03,P04,
        P05,P06,P07,P08,P09,
        P10,P11,P12,P13,P14,
        P15,P16,P17,P18,P19
    > > > &me, TIter it)
    {
    SEQAN_CHECKPOINT
        typedef HardwiredShape<
            P00,P01,P02,P03,P04,
            P05,P06,P07,P08,P09,
            P10,P11,P12,P13,P14,
            P15,P16,P17,P18,P19 >                                TSpec;
        typedef GappedShape<TSpec>                            TShape;
        typedef typename Value< Shape<TValue, TShape> >::Type    THValue;

        me.hValue = (THValue)ordValue((TValue)*it);
        return me.hValue = _hashHardwiredShape(me.hValue, it, TValue(), TSpec());
    }

//____________________________________________________________________________

    template <typename TValue, typename TSpec, typename TIter>
    inline void
    hashInit(Shape<TValue, GappedShape<TSpec> > &, TIter &)
    {
    }

    template <typename TValue, typename TSpec, typename TIter>
    inline typename Value< Shape<TValue, GappedShape<TSpec> > >::Type
    hashNext(Shape<TValue, GappedShape<TSpec> > &me, TIter &it)
    {
    SEQAN_CHECKPOINT
        return hash(me, it);
    }


//____________________________________________________________________________
    template <typename TValue, typename TSpec, typename TShapeString>
    inline bool
    stringToShape(
        Shape<TValue, GappedShape<TSpec> > &me,
        TShapeString const &bitmap)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TShapeString const>::Type        TIter;
        typedef typename Iterator<String<int> >::Type            TShapeIter;

        unsigned oneCount = 0;
        TIter it = begin(bitmap, Standard());
        TIter itEnd = end(bitmap, Standard());
        for(; it != itEnd; ++it)
            if (*it == '1')
                ++oneCount;

        if ((me.weight = oneCount) == 0) {
            me.span = 0;
            return false;
        }

        resize(me.diffs, oneCount - 1);

        unsigned diff = 1;
        me.span = 1;
        TShapeIter itS = begin(me.diffs, Standard());

        it = begin(bitmap, Standard());

        // skip leading zeros
        for(; it != itEnd && *it == '0' ; ++it) ;

        if (*it != '1')
            return false;

        for(++it; it != itEnd; ++it)
        {
            if (*it == '1')
            {
                me.span += diff;
                *itS = diff;
                ++itS;
                diff = 0;
            }
            ++diff;
        }

        return true;
    }

//____________________________________________________________________________

    template <typename TShapeString, typename TValue, typename TSpec>
    inline void
    shapeToString(
        TShapeString &bitmap,
        Shape<TValue, GappedShape<TSpec> > const &me)
    {
    SEQAN_CHECKPOINT

        clear(bitmap);
        if (weight(me) == 0) return;

        reserve(bitmap, length(me));
        appendValue(bitmap, '1');
        for (unsigned i = 0; i < weight(me) - 1; ++i)
        {
            for (int j = 1; j < me.diffs[i]; ++j)
                appendValue(bitmap, '0');
            appendValue(bitmap, '1');
        }
    }

//____________________________________________________________________________
    template <typename TValue, typename TSpec>
    inline void
    reverse(Shape<TValue, GappedShape<TSpec> > &me)
    {
    SEQAN_CHECKPOINT
        reverse(me.diffs);
    }

}    // namespace seqan

#endif
