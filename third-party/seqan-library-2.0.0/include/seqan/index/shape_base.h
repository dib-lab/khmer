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

 //TODO(weese:) sync new documentation to old dddoc updates

#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

    template <unsigned q>
    struct UngappedShape {};
    typedef UngappedShape<0> SimpleShape;

    template <typename TSpec>
    struct GappedShape {};
    typedef GappedShape<Default> GenericShape;


/*!
 * @class Shape
 * @headerfile <seqan/index.h>
 * @brief Stores hash value and shape for an ungapped or gapped q-gram.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class Shape;
 *
 * @tparam TSpec The specializing type. Default: @link SimpleShape @endlink, for
 *               ungapped q-grams.
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 *
 * The @link FiniteOrderedAlphabetConcept#ValueSize @endlink of Shape is the ValueSize of TValue which is the
 * alphabet size.
 *
 * To get the span or the weight of a shape call @link Shape#length @endlink or @link Shape#weight @endlink.
 */
    template <typename TValue = Dna, typename TSpec = SimpleShape>
    class Shape;

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn Shape#Value
 * @brief Returns the value type for a shape.
 *
 * @signature Value<TShape>::Type;
 *
 * @tparam TShape The Shape to query for its value type.
 *
 * @return Type The value type of the shape.
 */

    template <typename TValue, typename TSpec>
    struct Value<Shape<TValue,TSpec> >
    {
        typedef __uint64 Type;
    };

/*!
 * @mfn Shape#Size
 * @brief Returns the size type for a shape.
 *
 * @signature Size<TShape>::Type;
 *
 * @tparam TShape The Shape to query for its size type.
 *
 * @return Type The size type of the shape.
 */

    template <typename TValue, typename TSpec>
    struct Size<Shape<TValue,TSpec> >
    {
        typedef unsigned long Type;
    };

/*!
 * @mfn Shape#LENGTH
 * @brief Returns the length (span) of a shape.
 *
 * @signature LENGTH<TShape>::VALUE;
 *
 * @tparam TShape The Shape to query for its length (span).
 *
 * @return VALUE The length (span) of the shape.
 */

    template <typename TValue, unsigned q>
    struct LENGTH< Shape<TValue, UngappedShape<q> > >
    {
        enum { VALUE = q };
    };

/*!
 * @mfn Shape#WEIGHT
 * @brief Returns the weight (number of 1's) of a shape.
 *
 * @signature WEIGHT<TShape>::VALUE;
 *
 * @tparam TShape The Shape to query for its weight (number of 1's).
 *
 * @return VALUE The weight (number of 1's) of the shape.
 */

    template <typename TValue, unsigned q>
    struct WEIGHT< Shape<TValue, UngappedShape<q> > >
    {
        enum { VALUE = q };
    };

/*!
 * @mfn Shape#ValueSize
 * @brief Returns the type to use for the value size.
 *
 * @signature ValueSize<TShape>::Type;
 *
 * @tparam TShape The Shape to query for value size type.
 *
 * @return Type Type to use for the value size.
 */

    template <typename TValue, typename TSpec>
    struct ValueSize< Shape<TValue, TSpec> >
    {
        typedef typename Value<Shape<TValue, TSpec> >::Type THashValue;
        static const THashValue VALUE = Power<
                        ValueSize<TValue>::VALUE,
                        WEIGHT< Shape<TValue, TSpec> >::VALUE >::VALUE;
    };

/*!
 * @mfn Shape#Host
 * @brief Returns the host (= value) type to use.
 *
 * @signature Host<TShape>::Type;
 *
 * @tparam TShape The Shape to query for host (= value) type.
 *
 * @return Type Type to use for the host (= value) size.
 */

    template <typename TValue, typename TSpec>
    struct Host<Shape<TValue,TSpec> >
    {
        typedef TValue Type;
    };


//////////////////////////////////////////////////////////////////////////////

/*!
 * @class SimpleShape
 * @extends Shape
 * @headerfile <seqan/index.h>
 * @brief A variable length ungapped shape (also called q-gram or k-mer).
 *
 * @signature template <typename TValue>
 *            class Shape<TValue, SimpleShape>;
 *
 * @tparam TValue The @link Value @endlink type of the string the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 *
 * @section Remarks
 *
 * A SimpleShape must be resized first to a valid length.  To do so, call @link
 * Shape#resize @endlink.
 *
 * @see UngappedShape
 */

    //////////////////////////////////////////////////////////////////////////////
    // ungapped shape with variable length
    //////////////////////////////////////////////////////////////////////////////

    template <typename TValue>
    class Shape<TValue, SimpleShape>
    {
    public:
//____________________________________________________________________________

        unsigned                    span;
        typename Value<Shape>::Type    hValue;
        typename Value<Shape>::Type    XValue;
        typename Value<Shape>::Type    leftFactor;
        typename Value<Shape>::Type    leftFactor2;
        TValue                        leftChar;
//____________________________________________________________________________

/*!
 * @fn SimpleShape::Shape
 *
 * @brief Constructor.
 *
 * @signature Shape::Shape();
 * @signature Shape::Shape(shape);
 * @signature Shape::Shape(q);
 *
 * @param[in] shape Other Shape object (copy constructor).
 * @param[in] q     Length of the ungapped q-gram (@link IntegerConcept @endlink).
 */
        Shape():
            span(0),
            hValue(0),
            XValue(0),
            leftFactor(0),
            leftFactor2(0),
            leftChar(0) {}

        Shape(unsigned _span):
            hValue(0),
            XValue(0),
            leftFactor(0),
            leftFactor2(0),
            leftChar(0)
        {
            resize(*this, _span);
        }

        template <unsigned q>
        Shape(Shape<TValue, UngappedShape<q> > const &other)
        {
            *this = other;
        }

//____________________________________________________________________________

        template <unsigned q>
        inline Shape &
        operator=(Shape<TValue, UngappedShape<q> > const &other)
        {
            span = other.span;
            hValue = other.hValue;
            XValue = other.XValue;
            leftFactor = other.leftFactor;
            leftFactor2 = other.leftFactor2;
            leftChar = other.leftChar;
            return *this;
        }
    };

    //////////////////////////////////////////////////////////////////////////////
    // ungapped shape with fixed length q
    //////////////////////////////////////////////////////////////////////////////

/*!
 * @class UngappedShape
 * @extends Shape
 * @headerfile <seqan/index.h>
 * @brief A fixed length ungapped shape (also called q-gram or k-mer).
 *
 * @signature template <typename TValue, unsigned Q>
 *            class Shape<TValue, UngappedShape<Q> >;
 *
 * @tparam Q      The length of the shape (@link IntegerConcept @endlink).
 * @tparam TValue The @link Value @endlink type of the sequence the shape is
 *                applied to (e.g. <tt>Dna</tt>).
 *
 * @see SimpleShape
 */

    template <typename TValue, unsigned q>
    class Shape<TValue, UngappedShape<q> >
    {
    public:
        typedef typename Value<Shape>::Type THashValue;
//____________________________________________________________________________

        static const unsigned span = q;
        static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, q - 1>::VALUE;
        static const THashValue leftFactor2 = (Power<ValueSize<TValue>::VALUE, q>::VALUE - 1) / (ValueSize<TValue>::VALUE - 1);
        // Sigma^(q-1) + Sigma^(q-2) + ... + Sigma + 1

        THashValue    hValue;        // current hash value
        THashValue    XValue;        // Sum_{i=0..q-1} (x_i + 1)
        TValue        leftChar;    // leftmost character
//____________________________________________________________________________

        SEQAN_HOST_DEVICE
        Shape():
            hValue(0),
            XValue(0),
            leftChar(0) {}
    };



//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Shape#value
 * @headerfile <seqan/index.h>
 * @brief Returns the current hash value of the Shape.
 *
 * @signature TValue value(shape);
 *
 * @param[in] shape The Shape to query for its value.
 *
 * @return TValue The hash value of the shape.
 */

    template <typename TValue, typename TSpec>
    inline typename Value< Shape<TValue, TSpec> >::Type
    value(Shape<TValue, TSpec> &me)
    {
        return me.hValue;
    }

    template <typename TValue, typename TSpec>
    inline typename Value< Shape<TValue, TSpec> >::Type
    value(Shape<TValue, TSpec> const &me)
    {
        return me.hValue;
    }


//____________________________________________________________________________
/*!
 * @fn Shape#length
 * @brief Returns the number of elements of the shape (span).
 *
 * @signature TSize length(shape);
 *
 * @param[in] shape Shape object for which the number of relevant positions is determined.
 *
 * @return TSize The number of elements of the shape (span) (Metafunction: @link Shape#Size @endlink).
 */
    template <typename TValue, typename TSpec>
    inline SEQAN_HOST_DEVICE
    typename Size< Shape<TValue, TSpec> >::Type
    length(Shape<TValue, TSpec> const &me)
    {
    SEQAN_CHECKPOINT
        return me.span;
    }

//____________________________________________________________________________

/*!
 * @fn Shape#weight
 * @brief Number of relevant positions in a shape.
 *
 * @signature TSize weight(shape);
 *
 * @param[in] shape Shape object for which the number of relevant positions is determined.
 *
 * @return TSize Number of relevant positions  (Metafunction: @link Shape#Size @endlink).
 *
 * For ungapped shapes the return value is the result of the @link Shape#length
 * @endlink function.  For gapped shapes this is the number of '1's.
 */
    template <typename TValue, typename TSpec>
    inline SEQAN_HOST_DEVICE
    typename Size< Shape<TValue, TSpec> >::Type
    weight(Shape<TValue, TSpec> const &me)
    {
    SEQAN_CHECKPOINT
        return length(me);
    }

//____________________________________________________________________________

/*!
 * @fn Shape#resize
 * @brief Resize a shape to a specified span.
 *
 * @signature TSize resize(shape, length)
 *
 * @param[in,out] shape  Shape object for which the number of relevant positions is determined
 * @param[in]     length The new length (span) of the shape.
 *
 * @return TSize The new span of type  (Metafunction: @link Shape#Size @endlink).
 */
    template <typename TValue, typename TSize>
    inline typename Size< Shape<TValue, SimpleShape> >::Type
    resize(Shape<TValue, SimpleShape> & me, TSize new_length)
    {
    SEQAN_CHECKPOINT
        typedef typename Value< Shape<TValue, SimpleShape> >::Type    THValue;
        me.leftFactor = _intPow((THValue)ValueSize<TValue>::VALUE, new_length - 1);
        me.leftFactor2 = (_intPow((THValue)ValueSize<TValue>::VALUE, new_length) - 1) / (ValueSize<TValue>::VALUE - 1);
        return me.span = new_length;
    }

//____________________________________________________________________________

/*!
 * @fn Shape#hash
 * @brief Computes a (lower) hash value for a shape applied to a sequence.
 *
 * @signature TValue hash(shape, it[, charsLeft]);
 *
 * @param[in,out] shape     Shape to be used for hashing. Types: @link Shape @endlink
 * @param[in]     it        Sequence iterator pointing to the first character of the shape.
 * @param[in]     charsLeft The distance of <tt>it</tt> to the string end. If
 *                          <tt>charsLeft</tt> is smaller than the shape's span, the
 *                          hash value corresponds to the smallest shape beginning with
 *                          <tt>charsLeft</tt> characters.
 *
 * @return TValue Hash value of the shape  (Metafunction: @link Shape#Value @endlink).
 *
 * @see Shape#hashNext
 * @see Shape#hashUpper
 * @see Shape#hash2
 * @see Shape#unhash
 */
    template <typename TValue, typename TIter>
    inline typename Value< Shape<TValue, SimpleShape> >::Type
    hash(Shape<TValue, SimpleShape> &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, SimpleShape> >::Type    THValue;
        typedef typename Size< Shape<TValue, SimpleShape> >::Type    TSize;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.hValue = ordValue(me.leftChar = *it);
        for(TSize i = 1; i < me.span; ++i) {
            ++it;
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
        return me.hValue;
    }

/*!
 * @fn Shape#hashInit
 *
 * @brief Preprocessing step of a pure @link Shape#hashNext @endlink loop.
 *
 * Overlapping q-grams can efficiently be hashed by calling @link Shape#hash @endlink on the first text position and @link Shape#hashNext @endlink on succeeding, adjacent positions.
 * One drawback of this scenario is that for-loops cannot start with the first position directly and become more complicated.
 * As a remedy, @link Shape#hashInit @endlink was introduced which initializes the Shape to be used with @link Shape#hashNext @endlink on the first position directly.
 *
 * @signature void hashInit(shape, it);
 *
 * @param[in,out] shape Shape to be used for hasing.
 * @param[in]     it    The @link IteratorAssociatedTypesConcept iterator @endlink to use for initializing the shape.
 *
 * @section Example
 *
 * Two hash loop examples.
 * The first loop uses @link Shape#hash @endlink/@link Shape#hashNext @endlink while the second use @link Shape#hashInit @endlink/@link Shape#hashNext @endlink and can process all hashes within the loop.
 *
 * @include demos/index/shape_hash_init.cpp
 *
 * @code{.stdout}
 * 0    0    1    4    17    4    18    11    47    63    62    56
 * 0    0    1    4    17    4    18    11    47    63    62    56
 * @endcode
 */

    template <typename TValue, typename TIter>
    inline void
    hashInit(Shape<TValue, SimpleShape> &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, SimpleShape> >::Type    THValue;
        typedef typename Size< Shape<TValue, SimpleShape> >::Type    TSize;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.leftChar = 0;
        me.hValue = ordValue(*it);
        for(TSize i = 2; i < me.span; ++i) {
            ++it;
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
        }
    }

//____________________________________________________________________________
// fixed ungapped shapes

    // loop unrolling ...
    template <typename THValue, typename TValue, typename TIter>
    SEQAN_HOST_DEVICE inline THValue
    _hashFixedShape(THValue hash, TIter &, TValue const, UngappedShape<1> const) {
        return hash;
    }

    template <typename THValue, typename TValue, typename TIter, unsigned q>
    SEQAN_HOST_DEVICE inline THValue
    _hashFixedShape(THValue hash, TIter &it, TValue const, UngappedShape<q> const) {
        ++it;
        return _hashFixedShape(
            hash * ValueSize<TValue>::VALUE + ordValue((TValue)*it),
            it, TValue(), UngappedShape<q - 1>());
    }

    // ... for fixed ungapped shapes
    template <typename TValue, unsigned q, typename TIter>
    SEQAN_HOST_DEVICE inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
    hash(Shape<TValue, UngappedShape<q> > &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, UngappedShape<q> > >::Type    THValue;
        //typedef typename Size< Shape<TValue, UngappedShape<q> > >::Type     TSize;

        me.hValue = ordValue(me.leftChar = *it);
        return me.hValue = _hashFixedShape(me.hValue, it, TValue(), UngappedShape<q>());
    }

    template <typename TValue, unsigned q, typename TIter>
    inline typename Value< Shape<TValue, UngappedShape<q> > >::Type
    hashInit(Shape<TValue, UngappedShape<q> > &me, TIter it)
    {
        //typedef typename Value< Shape<TValue, UngappedShape<q> > >::Type    THValue;
        //typedef typename Size< Shape<TValue, UngappedShape<q> > >::Type    TSize;

        me.leftChar = 0;
        me.hValue = ordValue(*it);

        if (q > 1)
            me.hValue = _hashFixedShape(me.hValue, it, TValue(), UngappedShape<q-1>());

        return me.hValue;
    }

    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, TSpec> >::Type
    hash(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        TSize iEnd = me.span;
        if (iEnd > charsLeft) iEnd = charsLeft;

        TSize i = 0;
        if (iEnd > 0) {
            me.hValue = ordValue(me.leftChar = *it);
            for(i = 1; i < iEnd; ++i) {
                ++it;
                me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
            }
        } else
            return me.hValue = 0;

        // fill shape with zeros
        for(; i < (TSize)me.span; ++i)
            me.hValue *= ValueSize<TValue>::VALUE;
        return me.hValue;
    }

//____________________________________________________________________________
// Tuple -> fixed ungapped shapes

    template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TPack>
    inline THValue
    _hashTuple2FixedShape(
        THValue const,
        Tuple<TTValue, SIZE, TPack> const &tuple,
        TValue const,
        UngappedShape<1> const)
    {
        return ordValue(tuple[0]);
    }

    template <typename THValue, typename TValue, typename TTValue, unsigned SIZE, typename TPack, unsigned q>
    inline THValue
    _hashTuple2FixedShape(
        THValue const,
        Tuple<TTValue, SIZE, TPack> const &tuple,
        TValue const,
        UngappedShape<q> const)
    {
        return _hashTuple2FixedShape(THValue(), tuple, TValue(), UngappedShape<q - 1>())
            * ValueSize<TValue>::VALUE + ordValue(tuple[q-1]);
    }

    // ... for fixed ungapped shapes
    template <
        typename TValue,
        typename TTValue,
        unsigned SIZE,
        unsigned q>
    typename Value< Shape<TValue, UngappedShape<q> > >::Type
    hash(
        Shape<TValue, UngappedShape<q> > &me,
        Tuple<TTValue, SIZE, BitPacked<> > /*const &*/tuple)
    {
    SEQAN_CHECKPOINT
        if (ValueSize<TValue>::VALUE == (1 << BitsPerValue<TTValue>::VALUE))
            if (q == SIZE)
                return tuple.i;
            else
                return tuple >> (q - SIZE);
        else
            return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), UngappedShape<q>());
    }

    template <
        typename TValue,
        typename TTValue,
        unsigned SIZE,
        typename TPack,
        unsigned q>
    typename Value< Shape<TValue, UngappedShape<q> > >::Type
    hash(
        Shape<TValue, UngappedShape<q> > &me,
        Tuple<TTValue, SIZE, TPack> /*const &*/tuple)
    {
    SEQAN_CHECKPOINT
        return me.hValue = _hashTuple2FixedShape(me.hValue, tuple, TValue(), UngappedShape<q>());
    }

//____________________________________________________________________________

/*!
 * @fn Shape#hashUpper
 * @brief Computes an upper hash value for a shape applied to a sequence.
 *
 * @signature TValue hashUpper(shape, it, charsLeft);
 *
 * @param[in,out] shape     Shape to be used for hashing. Types: @link Shape @endlink
 * @param[in]     it        Sequence iterator pointing to the first character of the shape.
 * @param[in]     charsLeft The distance of <tt>it</tt> to the string end.
 *
 * @return TValue Upper hash value of the shape. The hash value corresponds to
 *                the maximal @link Shape#hash @endlink value of a shape beginning
 *                with <tt>min(charsLeft,length(shape))</tt> characters + 1  (Metafunction:
 *                @link Shape#Value @endlink).
 *
 * This function in conjunction with @link Shape#hash @endlink is useful to search a q-gram index for p-grams with
 * p &lt; q.
 *
 * @see Shape#hash
 */

    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, TSpec> >::Type
    hashUpper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        TSize iEnd = me.span;
        if (iEnd > charsLeft) iEnd = charsLeft;

        TSize i = 0;
        if (iEnd > 0) {
            me.hValue = ordValue(me.leftChar = *it);
            for(i = 1; i < iEnd; ++i) {
                ++it;
                me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*it);
            }
            ++me.hValue;
        } else
            me.hValue = 1;

        // fill shape with zeros
        for(; i < (TSize)me.span; ++i)
            me.hValue *= ValueSize<TValue>::VALUE;
        return me.hValue;
    }

//____________________________________________________________________________

/*!
 * @fn Shape#hashNext
 * @headerfile <seqan/index.h>
 * @brief Computes the hash value for the adjacent shape.
 *
 * @signature TValue hashNext(shape, it);
 *
 * @param[in,out] shape Shape to be used for hashing. Types: @link Shape @endlink
 * @param[in]     it    Sequence iterator pointing to the first character of the adjacent shape.
 *
 * @return TValue Hash value of the q-gram (Metafunction: @link Shape#Value @endlink).
 *
 * @link Shape#hash @endlink has to be called before.
 *
 * @see Shape#hash
 */
    template <typename TValue, typename TSpec, typename TIter>
    inline typename Value< Shape<TValue, TSpec> >::Type
    hashNext(Shape<TValue, TSpec> &me, TIter const &it)
    {
    SEQAN_CHECKPOINT
        // remove first, shift left, and add next character
        typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;
        typedef typename Size< Shape<TValue, TSpec> >::Type        TSize;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.hValue =
            (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE
            + ordValue((TValue)*(it + ((TSize)me.span - 1)));
        me.leftChar = *it;
        return me.hValue;
    }

//____________________________________________________________________________

/*!
 * @fn Shape#hash2
 * @brief Computes an unique hash value of a shape applied to a sequence, even if the sequence is shorter than
 *        the shape span.
 *
 * @signature TValue hash2(shape, it, charsLeft);
 *
 * @param[in,out] shape     Shape to be used for hashing. Types: @link Shape @endlink
 * @param[in]     it        Sequence iterator pointing to the first character of the shape.
 * @param[in]     charsLeft The distance of <tt>it</tt> to the string end.
 *
 * @return TValue Hash value of the shape (Metafunction: @link Shape#Value @endlink).
 *
 * @see Shape#hash2Next
 * @see Shape#hash2Upper
 * @see Shape#unhash
 * @see Shape#hash
 */

    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, TSpec> >::Type
    hash2(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
    {
        //typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        TSize iEnd = me.span;
        if (iEnd > charsLeft) iEnd = charsLeft;

        TSize i = 0;
        if (iEnd > 0) {
            me.hValue = me.XValue = ordValue(me.leftChar = *it);
            for(i = 1; i < iEnd; ++i) {
                ++it;
                // update sum of x_i
                me.XValue += ordValue((TValue)*it);
                // shift hash
                me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
            }
        } else
            return me.hValue = me.XValue = 0;

        // fill shape with zeros
        for(; i < (TSize)me.span; ++i)
            me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
        return me.hValue += iEnd;
    }

/*!
 * @fn Shape#hash2Upper
 * @brief Computes an upper unique hash value of a shape applied to a sequence,
 *        even if the sequence is shorter than the shape span.
 *
 * @signature TValue hash2Upper(shape, it, charsLeft);
 *
 * @param[in] shape     Shape to be used for hashing. Types: @link Shape @endlink
 * @param[in] it        Sequence iterator pointing to the first character of the shape.
 * @param[in] charsLeft The distance of <tt>it</tt> to the string end.
 *
 * @return TValue Upper hash value of the shape. The hash value corresponds to
 *                the maximal @link Shape#hash2 @endlink value of a shape beginning
 *                with the <tt>min(charsLeft,length(shape))</tt> characters + 1
 *                (Metafunction: @link Shape#Value @endlink).
 *
 * This function in conjunction with @link Shape#hash2 @endlink is useful to search a
 * q-gram index for p-grams with <i>p &lt; q</i>.
 *
 * @see Shape#hash2
 */
    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, TSpec> >::Type
    hash2Upper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
    {
    SEQAN_CHECKPOINT
        typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        TSize iEnd = me.span;
        if (iEnd > charsLeft) iEnd = charsLeft;

        THValue hValue, XValue;
        TSize i = 0;
        if (iEnd > 0) {
            hValue = XValue = ordValue((TValue)*it);
            for(i = 1; i < iEnd; ++i) {
                ++it;
                // update sum of x_i
                XValue += ordValue((TValue)*it);
                // shift hash
                hValue = hValue * ValueSize<TValue>::VALUE + XValue;
            }
        } else
            hValue = XValue = 0;

        if (charsLeft <= me.span) {
            ++XValue;
            ++hValue;
        }

        // fill shape with zeros
        for(; i < (TSize)me.span; ++i)
            hValue = hValue * ValueSize<TValue>::VALUE + XValue;
        return hValue += iEnd;
    }

//____________________________________________________________________________

/*!
 * @fn Shape#hash2Next
 * @headerfile <seqan/index.h>
 * @brief Computes a unique hash value for the adjacent shape, even if it is shorter than q.
 *
 * @signature TValue hash2Next(shape, it);
 *
 * @param[in,out] shape Shape to be used for hashing the q-gram. Types: @link Shape @endlink
 * @param[in,out] it    Sequence iterator pointing to the first character of the adjacent shape.
 *
 * @return TValue Hash value of the shape (Metafunction: @link Shape#Value @endlink).
 *
 * @link Shape#hash @endlink has to be called before with <tt>shape</tt> on the left adjacent q-gram.
 *
 * @see Shape#hash2
 */

    template <typename TValue, typename TSpec, typename TIter, typename TSize>
    inline typename Value< Shape<TValue, TSpec> >::Type
    hash2Next(Shape<TValue, TSpec> &me, TIter &it, TSize charsLeft)
    {
    SEQAN_CHECKPOINT
        // remove first, shift left, and add next character
        typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        if (charsLeft >= me.span) {
            // update sum of x_i
            me.XValue = me.XValue + ordValue((TValue)*(it + me.span - 1)) - ordValue(me.leftChar);
            // shift hash
            me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
                        - me.span * (ValueSize<TValue>::VALUE - 1);
        } else {
            // update sum of x_i
            me.XValue -= ordValue(me.leftChar);
            // shift hash
            me.hValue = (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
                        - charsLeft * (ValueSize<TValue>::VALUE - 1) - ValueSize<TValue>::VALUE;
        }

        me.leftChar = *it;
        return me.hValue;
    }

/*!
 * @fn Shape#unhash
 * @headerfile <seqan/index.h>
 * @brief Inverse of the @link Shape#hash @endlink function; for ungapped shapes.
 *
 * @signature void unhash(result, hash, q);
 *
 * @param[out] result @link String @endlink to write the result to. Types: @link String @endlink.
 * @param[in]  hash   The hash value previously computed with @link Shape#hash @endlink.
 * @param[in]  q      The <tt>q</tt>-gram length. Types: <tt>unsigned</tt>
 *
 * @see Shape#hash
 * @see Shape#hash2
 */
    template <typename TString, typename THash>
    inline void unhash(TString &result, THash hash, unsigned q)
    {
    SEQAN_CHECKPOINT
        typedef typename Value<TString>::Type    TValue;

        resize(result, q);
        for (unsigned i = q; i > 0; )
        {
            result[--i] = (TValue)(hash % ValueSize<TValue>::VALUE);
            hash /= ValueSize<TValue>::VALUE;
        }
    }

//____________________________________________________________________________

/*!
 * @fn Shape#stringToShape
 * @brief Takes a shape given as a string of '1' (relevant position) and '0'
 *        (irrelevant position) and converts it into a Shape object.
 *
 * @signature bool stringToShape(shape, bitmap);
 *
 * @param[in,out] shape Shape object that is manipulated.
 * @param[in]     bitmap A character string of '1' and '0' representing relevant and irrelevant positions (blanks)
 *                respectively.  This string must begin with a '1'.  Trailing '0's are ignored.  If <tt>shape</tt>
 *                is a @link SimpleShape @endlink at most one contiguous sequences of <tt>1</tt>s is allowed.  If
 *                <tt>shape</tt> is a @link OneGappedShape @endlink at most two contiguous sequences of '1's are
 *                allowed (@link String @endlink of <tt>char</tt>).
 *
 * @return bool <tt>true</tt> if the conversion was successful.
 *
 * @see Shape#shapeToString
 * @see reverse
 */
    template <typename TValue, typename TShapeString>
    inline bool
    stringToShape(
        Shape<TValue, SimpleShape> &me,
        TShapeString const &bitmap)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TShapeString const>::Type        TIter;
        typedef typename Size<TShapeString const>::Type            TSize;

        TIter it = begin(bitmap, Standard());
        TIter itEnd = end(bitmap, Standard());

        TSize ones = 0;
        for(; it != itEnd && *it == '0' ; ++it) ;
        for(; it != itEnd && *it == '1' ; ++it)    ++ones;
        for(; it != itEnd && *it == '0' ; ++it) ;

        resize(me, ones);

        return it == itEnd;
    }

//____________________________________________________________________________

/*!
 * @fn Shape#shapeToString
 * @brief Converts a given shape into a sequence of '1' (relevant position) and '0' (irrelevant position).
 *
 * @signature void shapeToString(bitmap, shape);
 *
 * @param[in,out] bitmap The resulting sequence object. Types: @link String @endlink
 * @param[in]     shape Shape object. Types: @link Shape @endlink
 *
 * @see Shape#stringToShape
 */

    template <typename TShapeString, typename TValue, unsigned q>
    inline void
    shapeToString(
        TShapeString &bitmap,
        Shape<TValue, UngappedShape<q> > const &me)
    {
    SEQAN_CHECKPOINT

        clear(bitmap);
        resize(bitmap, length(me), '1');
    }

//____________________________________________________________________________

    template <typename TValue, typename TSpec>
    inline void
    reverse(Shape<TValue, TSpec> &)
    {
    }

}    // namespace seqan

#endif
