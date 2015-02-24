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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================
// Definition of generic and fixed-at-compile-time CyclicShapes
// ==========================================================================

#ifndef SEQAN_HEADER_CYCLIC_SHAPE_H
#define SEQAN_HEADER_CYCLIC_SHAPE_H

#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// needed for GappedShape<HardwiredShape<...> >
template<typename TSpec> struct GappedShape;

// GenericShape defined in shape_gapped.h via typedef
typedef GappedShape<Default> GenericShape;

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Class CyclicShape
// --------------------------------------------------------------------------

/*!
 * @class CyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief A pattern of zeros and ones to mark "don't care"-positions in a text.
 *
 * @signature template <typename TShapeSpec>
 *            class CyclicShape<TShapeSpec>;
 *
 * @tparam TShapeSpec The specializing type. Default is @link GenericCyclicShape
 *         GenericShape @endlink, another option is @link FixedCyclicShape
 *         FixedShape @endlink.
 * @tparam TSize Size type of the CyclicShape
 *
 * CyclicShapes store a pattern like <tt>11010100110</tt> that mark
 * care (match) positions and don't-care positions.
 * Unlike @link Shape @endlink, the CyclicShape does not perform hashing of q-grams.
 * It is instead useful to modify a text in such a way that zero
 * position are ignored. The pattern is applied repeatedly on the whole
 * text, as the example shows.
 *
 * Note that CyclicShapes can start and end with zero characters, which
 * is not allowed in @link Shape @endlink.
 *
 * @section Examples
 * @include demos/cyclic_shape.cpp
 *
 * The output is as follows:
 *
 * @include demos/cyclic_shape.cpp.stdout
 *
 * @see ModCyclicShapeModifiedString
 * @see ModCyclicShapeModifiedIterator
 */
template<typename TCyclicSpec = GenericShape>
class CyclicShape;                           // member: diffs, loffset, span

template<unsigned L, typename THardWiredShape, unsigned R>
struct FixedShape;

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

/*!
 * @mfn CyclicShape#Size
 * @headerfile <seqan/modifier.h>
 *
 * @brief Size type for parameters used in CyclicShape.
 *
 * @signature Size<CyclicShape<TSpec> >::Type;
 *
 * @tparam TSpec The CyclicShape specialization.
 *
 * @return TReturn Currently the return type <tt>unsigned char</tt>.
 *
 * @section Remarks
 */
template<typename TSpec>
struct Size<CyclicShape<TSpec> >
{
    typedef unsigned char Type;
};

// --------------------------------------------------------------------------
// Class GenericCyclicShape
// --------------------------------------------------------------------------

/*!
 * @class GenericCyclicShape Generic CyclicShape
 * @extends CyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief CyclicShape initiated at runtime.
 *
 * @tparam TSize Size type of the CyclicShape
 *
 * @signature template <>
 *            class CyclicShape<GenericShape>;
 *
 * Generic Shapes can be easily created using <tt>stringToCyclicShape</tt>,
 * but their usage comes with a significant loss in efficiency.
 * For longer sequences, prefer the @link FixedCyclicShape @endlink.
 *
 * @section Examples
 * @include demos/cyclic_shape.cpp
 *
 * The output is as follows:
 *
 * @include demos/cyclic_shape.cpp.stdout
 *
 * @see ModCyclicShapeModifiedString
 * @see ModCyclicShapeModifiedIterator
 */
template<>
class CyclicShape<GenericShape>
{
public:
    /*!
     * @var String<TSize> GenericCyclicShape::diffs;
     * @brief Distances between care positions.
     *
     * <tt>diffs</tt> has <i>weight</i> many non-zero entries describing the size
     * of the gap after a care-position. The last entry holds the cyclic distance from
     * the last "1" to the first one.
     */
    String<Size<CyclicShape>::Type> diffs;

    /*!
     * @var TSize CyclicShape::span;
     * @brief span of the CyclicShape.
     * @var TSize CyclicShape::loffset;
     * @brief left offset (number of leading zeros) of the CyclicShape.
     */
    Size<CyclicShape>::Type loffset, span;

    /*!
     * @fn GenericCyclicShape::CyclicShape
     * @brief The constructor
     * @signature CyclicShape();
     * @signature CyclicShape(shape);
     * @param[in] shape A GenericCyclicShape to be copied.
     *
     * The default constructor generates the pattern "1".
     * The copy constructor copies the GenericCyclicShape shape.
     *
     * @snippet demos/cyclic_shape_snippets.cpp Define GenericCyclicShape
     *
     * @Remarks
     *
     * If you want to generate a GenericCyclicShape from a FixedCyclicShape,
     * use @link CyclicShape#cyclicShapeToString @endlink and @link GenericCyclicShape#stringToCyclicShape @endlink
     */
    CyclicShape() :
        loffset(0), span(1)
    {
        // Pattern "1"
        appendValue(diffs, 1);
    }

    CyclicShape(CyclicShape const & other) :
        diffs(other.diffs), loffset(other.loffset), span(other.span)
    {}
};

// --------------------------------------------------------------------------
// Class FixedCyclicShape
// --------------------------------------------------------------------------

/*!
 * @class FixedCyclicShape Fixed CyclicShape
 * @extends CyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief CyclicShape defined at compile time.
 *
 * @signature template <unsigned L, typename THardwiredShape, unsigned R>
 *            class CycledShape<FixedShape<L, GappedShape<THardwiredShape>, R> >;
 *
 * @tparam L Left offset. Number of leading zeros of the Shape.
 * @tparam R Right offset. Number of trailing zeros of the Shape.
 * @tparam THardwiredShape A specialization of @link HardwiredShape @endlink
 * @tparam TSize Size type of the CyclicShape
 *
 * Fixed CylcicShapes contain their information at compile time, so in
 * most cases no object must be created.
 * The notation is similar to the one of @link HardwiredShape @endlink: Template
 * paramters code for the distance between ones. Additionally you have to specify
 * the number of leading and trailing zeros.
 *
 * The notation is chosen in such a way that predefined Shapes like
 * PatternHunter can be plugged into a CyclicShape. Like
 * HardwiredShapes, Fixed CyclicShapes are limited to a weight of 21.
 *
 * See @link CyclicShape @endlink for an example on how to use a CyclicShape.
 *
 * @see ModCyclicShapeModifiedString
 * @see ModCyclicShapeModifiedIterator
 */
template<unsigned LeftOffset, unsigned RightOffset,
         int P00, int P01, int P02, int P03, int P04,
         int P05, int P06, int P07, int P08, int P09,
         int P10, int P11, int P12, int P13, int P14,
         int P15, int P16, int P17, int P18, int P19>
class CyclicShape<FixedShape<LeftOffset, GappedShape<HardwiredShape<
        P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19>
        >, RightOffset> >
{
public:
    typedef HardwiredShape<
            P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19>
            THardWiredShape;

    enum {loffset = LeftOffset};
    enum {span = LeftOffset + LENGTH<THardWiredShape>::VALUE + RightOffset};

    /*!
     * @var TSize[] FixedCyclicShape::diffs;
     * @brief Distances between care positions.
     *
     * <tt>static const TSize diffs[]</tt> has <i>weight</i> many non-zero entries
     * in the beginning, the remaining ones are zero. An additional entry at position
     * <i>weight1</i> holds the cyclic distance from the last "1" to the first one.
     * During the iteration with a @link ModCyclicShapeModifiedIterator @endlink
     * an index position keeps track of the position inside <tt>diffs</tt>.
     */
    static const typename Size<CyclicShape>::Type diffs[];

    /*!
     * @fn FixedCyclicShape::CyclicShape
     * @brief The constructor
     *
     * @signature CyclicShape();
     * @signature CyclicShape(shape);
     *
     * @param[in] shape A FixedCyclicShape of this type
     *
     * This constructor does not do anything, the FixedCyclicShape is defined by its type alone.
     * See the example on how to create a CyclicShape:
     *
     * @snippet demos/cyclic_shape_snippets.cpp Define FixedCyclicShape
     */
    CyclicShape()
    {}

    CyclicShape(CyclicShape const &)
    {}
};

// --------------------------------------------------------------------------
// static const int CyclicShape<FixedShape<L, TSpec, R> >::diffs[]
//
// ugly work around to derive the static const array diffs[] from
// the static const array HardwiredShape::DIFFS. In diffs, the first
// zero-entry of DIFFS had to be changed.
// Note: This requires that the non-used entries of
// HardwiredShape::DIFFS are zero.
// --------------------------------------------------------------------------

template<unsigned entry, typename TFixedS>
struct _tmpArrEntry;

// normal positions in diff[]
template<unsigned entry, unsigned L, unsigned R,
         int P00, int P01, int P02, int P03, int P04,
         int P05, int P06, int P07, int P08, int P09,
         int P10, int P11, int P12, int P13, int P14,
         int P15, int P16, int P17, int P18, int P19>
struct _tmpArrEntry<entry, FixedShape<L, GappedShape<HardwiredShape
        <P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19> >, R> >
{
    enum {VALUE = entry};
};

// zero positions in diff[]
template<unsigned L, unsigned R,
         int P00, int P01, int P02, int P03, int P04,
         int P05, int P06, int P07, int P08, int P09,
         int P10, int P11, int P12, int P13, int P14,
         int P15, int P16, int P17, int P18, int P19>
struct _tmpArrEntry<0, FixedShape<L, GappedShape<HardwiredShape
        <P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19> >, R> >
{
    enum {VALUE = L + R + 1};
};

#define ENTRY(X) \
    _tmpArrEntry<X, FixedShape<L, GappedShape<HardwiredShape \
    <P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, \
    P10, P11, P12, P13, P14, P15, P16, P17, P18, P19> >, R> >::VALUE

template<unsigned L, unsigned R,
        int P00, int P01, int P02, int P03, int P04,
        int P05, int P06, int P07, int P08, int P09,
        int P10, int P11, int P12, int P13, int P14,
        int P15, int P16, int P17, int P18, int P19>
const typename Size<CyclicShape<FixedShape<L, GappedShape<HardwiredShape
        <P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19>
        >, R> > >::Type
CyclicShape<FixedShape<L, GappedShape<HardwiredShape
        <P00, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19>
        >, R> >::diffs[] =
{
    ENTRY(P00), ENTRY(P01), ENTRY(P02), ENTRY(P03), ENTRY(P04),
    ENTRY(P05), ENTRY(P06), ENTRY(P07), ENTRY(P08), ENTRY(P09),
    ENTRY(P10), ENTRY(P11), ENTRY(P12), ENTRY(P13), ENTRY(P14),
    ENTRY(P15), ENTRY(P16), ENTRY(P17), ENTRY(P18), ENTRY(P19),
    ENTRY(0)
};

#undef ENTRY



// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction WEIGHT
// --------------------------------------------------------------------------

/*!
 * @mfn FixedCyclicShape#WEIGHT
 * @headerfile <seqan/modifier.h>
 *
 * @brief Weight (number of care-positions) of Fixed CyclicShapes
 *
 * @signature WEIGHT<CyclicShape<FixedShape<L, TInnerShape, R> > >::VALUE;
 *
 * @tparam L left offset
 * @tparam R right offset
 * @tparam TInnerShape specialization of a hardwired GappedShape
 *
 * @return Weight of a fixed CyclicShape (enum).
 *
 * @see CyclicShape#weight
 */

template<unsigned L, typename THardwiredShape, unsigned R>
struct WEIGHT<CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > >
{
    enum {VALUE = WEIGHT<THardwiredShape>::VALUE};
};


// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function weight
// --------------------------------------------------------------------------

/*!
 * @fn CyclicShape#weight
 * @headerfile <seqan/modifier.h>
 *
 * @brief Return the weight of a CyclicShape
 *
 * @signature TSize weight(const & cyclicShape);
 *
 * @tparam TSpec Specialisation of the CyclicShape.
 * @tparam TSize Size type of CyclicShape
 * @param[in] cyclicShape CyclicShape object
 * @return weight of the CyclicShape (number of care positions)
 */
template<typename TSpec>
inline typename Size<CyclicShape<TSpec> >::Type
weight(CyclicShape<TSpec> const & s)
{
    return static_cast<Size<CyclicShape<> >::Type>(length(s.diffs));
}

template<unsigned L, typename THardwiredShape, unsigned R>
inline typename Size<CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > >::Type
weight(CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > const &)
{
    typedef CyclicShape<FixedShape<L, GappedShape<THardwiredShape>, R> > TCyclicShape;
    return WEIGHT<TCyclicShape>::VALUE;
}

// --------------------------------------------------------------------------
// Function stringToCyclicShape()
// --------------------------------------------------------------------------

/*!
 * @fn GenericCyclicShape#stringToCyclicShape
 * @headerfile <seqan/modifier.h>
 *
 * @brief Converts a 0/1 string to a Generic CyclicShape.
 *
 * @signature void stringToCyclicShape(shape, bitmap);
 *
 * @tparam TString A string type, e.g. CharString.
 * @param[out] shape Generic CyclicShape
 * @param[in] bitmap 0/1 string of type TString. CyclicShapes may start and end with zeros,
 *            but must contain at least one 1.
 *
 */
template<typename TString>
inline bool
stringToCyclicShape(CyclicShape<GenericShape> & shape, TString const & bitmap)
{
    typename Iterator<TString const>::Type it    = begin(bitmap);
    typename Iterator<TString const>::Type itBeg = begin(bitmap);
    typename Iterator<TString const>::Type itEnd = end(bitmap);

    Size<CyclicShape<GenericShape> >::Type countOnes = 0;
    for (; it != itEnd; ++it)
    {
        if (*it == '1')
            ++countOnes;
    }
    SEQAN_ASSERT_GT(countOnes, 0u);

    resize(shape.diffs, countOnes);

    typename Iterator<String<Size<CyclicShape<GenericShape> >::Type> >::Type diffIter = begin(shape.diffs);

    // go to first 1
    for (it = itBeg; *it != '1'; ++it)
        continue;
    shape.loffset = (Size<CyclicShape<GenericShape> >::Type)(it - itBeg);

    countOnes = 1; // now used as the distance between 1s
    for (++it; it != itEnd; ++it)
    {
        if (*it == '1')
        {
            *diffIter = countOnes;
            ++diffIter;
            countOnes = 1;
        }
        else
            ++countOnes;
    }

    *diffIter = countOnes + shape.loffset;
    shape.span = (typename Size<CyclicShape<GenericShape> >::Type)(it - itBeg);
    return true;
}

// --------------------------------------------------------------------------
// Function cyclicShapeToString()
// --------------------------------------------------------------------------

/*!
 * @fn CyclicShape#cyclicShapeToString
 *
 * @brief Print a given cyclic shape as a sequence of '1' (care position) and
 *        '0' (don't-care position).
 *
 * @signature void cyclicShapeToString(bitmap, cyclicShape);
 *
 * @param[in] cyclicShape CyclicShape object. Types: @link CyclicShape @endlink
 * @param[out] bitmap The resulting sequence object. Type: @link String @endlink,
 *             e.g. CharString
 *
 * @see GenericCyclicShape#stringToCyclicShape
 */
template<typename TShapeString, typename TSpec>
inline void
cyclicShapeToString(
    TShapeString & bitmap,
    CyclicShape<TSpec> const & me)
{
    clear(bitmap);
    if (weight(me) == 0)
        return;

    typename Size<CyclicShape<TSpec> >::Type i, j;
    for (i = 0; i < me.loffset; ++i)
        appendValue(bitmap, '0');
    for (i = 0; i < weight(me) - 1; ++i)
    {
        appendValue(bitmap, '1');
        for (j = 1; j < static_cast<typename Size<CyclicShape<TSpec> >::Type>(me.diffs[i]); ++j)
            appendValue(bitmap, '0');
    }
    appendValue(bitmap, '1');
    for (i = length(bitmap); i < me.span; ++i)
        appendValue(bitmap, '0');
}

// --------------------------------------------------------------------------
// Function carePositions()
// --------------------------------------------------------------------------

/*!
 * @fn CyclicShape#carePositions
 *
 * @brief Determine the indices of care positions in the range <i>[0,span)</i>
 *
 * @signature void carePositions(positions, cyclicShape);
 *
 * @param[in] cyclicShape CyclicShape object.
 * @param[out] positions The resulting @link String @endlink to store the care-positions in. Type TString.
 *
 * @tparam TString String or array. Value type should be an integral size type.
 * @tparam TSpec Specialization of CyclicShape.
 *
 * This function can be used to convert <tt>CyclicShape.diffs</tt>, which stores
 * the distances between care positions, to a <tt>positions</tt> string directly
 * containing the care positions. See the example:
 *
 * @snippet demos/cyclic_shape_snippets.cpp CyclicShape Care Positions
 */

template<typename TString, typename TSpec>
inline void
carePositions(TString & output, CyclicShape<TSpec> const & shape)
{
    typedef typename Size<CyclicShape<TSpec> >::Type TPos;

    resize(output, weight(shape));
    TPos val = shape.loffset;
    output[0] = shape.loffset;

    for (TPos i = 1; i < weight(shape); ++i)
    {
        val += shape.diffs[i - 1];
        output[i] = val;
    }
}

// --------------------------------------------------------------------------
// Function cyclicShapeToSuffixLengths()
// --------------------------------------------------------------------------

/*!
 * @fn CyclicShape#cyclicShapeToSuffixLengths
 *
 * @brief Calculates the number of characters of modified Strings
 *        shorter than the Shape's span.
 *
 * @signature void cyclicShapeToSuffixLengths(TString suffLengths, cyclicShape);
 *
 * @param[in] cyclicShape CyclicShape object.
 * @param[in,out] suffLengths String to be filled. Must be resized beforehands.
 *                Fixed length arrays also work.
 *
 * @tparam TString String type or array with an integral value type, i.e. unsigned or int.
 *
 * Given a CyclicShape, this function calculates a how many characters a String of length
 * <i>x</i> contains after the CyclicShape is applied to it. This is done for all <tt>0 <= x < span</tt>.
 * <tt>suffLengths</tt> therefore must be resized to the shape's span beforehands.
 * The resizing is not done in this function so that it can be applied to arrays, too.
 */

template<typename TString,
         typename TShape>
inline void cyclicShapeToSuffixLengths(TString & suffLengths, TShape const & shape)
{
    typedef typename Size<TString>::Type TSize;

    TSize w = weight(shape);
    TSize s = shape.span;
    TSize o = shape.loffset;

    SEQAN_ASSERT_GEQ(s, 1u);
    //SEQAN_ASSERT_EQ(length(suffLengths), s); // Disabled to support arrays too.

    // build cummulative some of distances first
    String<TSize> sums;
    resize(sums, w);
    sums[0] = o;
    for (unsigned i = 0; i < w - 1; ++i)
        sums[i + 1] = sums[i] + shape.diffs[i];

    // write suffixLengths
    unsigned sumPos = 0;
    for (unsigned i = 0; i < s; ++i)
        suffLengths[i] = (sumPos < w && i == sums[sumPos]) ? sumPos++ : sumPos;
}

} //namespace


#endif
