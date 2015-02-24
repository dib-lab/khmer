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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Switch to Host interface.

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec>
struct GapsIterator;

/*!
 * @defgroup GapsSpecTag Gaps Specialization Tags
 * @brief Tags for specializing the Gaps class.
 */

/*!
 * @tag GapsSpecTag#ArrayGaps
 * @headerfile <seqan/align.h>
 * @brief Tag for the Array Gaps specialization.
 *
 * @signature struct ArrayGaps_;
 * @signature typedef Tag<ArrayGaps_> ArrayGaps;
 */

struct ArrayGaps_;
typedef Tag<ArrayGaps_> ArrayGaps;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Gaps
// ----------------------------------------------------------------------------

/*!
 * @class Gaps
 * @implements ContainerConcept
 * @headerfile <seqan/align.h>
 * @brief Store the gapped version of a sequence.
 *
 * @signature template <typename TSequence, typename TSpec>
 *            class Gaps;
 *
 * @tparam TSequence The type of the underlying sequence.
 * @tparam TSpec     Tag for specialization.
 *
 * Gaps wrap a @link ContainerConcept Sequence @endlink and allows to (1) insert gaps into the sequence and (2) select
 * an infix of the gapped sequence (clipping).  The gaps are not inserted into the underlying sequence (source) but
 * stored separately.  Using the clipping is optional and meant for selecting parts of the alignment as a part of the
 * result of a local alignment algorithm.
 *
 * <img src="gaps_illustration.png" title="Illustration of Gaps object and positions with clipping." />
 *
 * In the figure above, the source sequence has seven characters, the gapped sequence has four gaps and thus consists
 * of eleven characters.  The gapped sequence is clipped to start at position 0 in the gapped sequence and to end at
 * position 8 in the gapped sequence (the positions given as half-open intervals <tt>[begin, end)</tt>).
 *
 * The figure shows the three coordinate systems that are used with Gaps objects.  The source position is the position
 * in the underlying sequence.  The unclipped view position is the position in the gapped sequence without gaps.  The
 * view position is the position in the gapped sequence but including the clipping: All (clipped) view positions have
 * the clipping begin position subtracted from them.
 *
 * @section Examples
 *
 * The following example shows the construction of the gaps object from the image above together with some calls to
 * <tt>toViewPosition</tt> and <tt>toSourcePosition</tt>.
 *
 * @include demos/align/gaps_example.cpp
 *
 * The output is as follows:
 *
 * @include demos/align/gaps_example.cpp.stdout
 */

template <typename TSequence, typename TSpec = ArrayGaps>
class Gaps;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
struct Value<Gaps<TSequence, TSpec> >
{
    typedef typename Value<TSequence>::Type           TAlphabet;
    typedef typename GappedValueType<TAlphabet>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Value<Gaps<TSequence, TSpec> const> : Value<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec, typename TIteratorSpec>
struct Iterator<Gaps<TSequence, TSpec>, TIteratorSpec>
{
    typedef Iter<Gaps<TSequence, TSpec>, GapsIterator<TSpec> > Type;
};

template <typename TSequence, typename TSpec, typename TIteratorSpec>
struct Iterator<Gaps<TSequence, TSpec> const, TIteratorSpec>
{
    typedef Iter<Gaps<TSequence, TSpec> const, GapsIterator<TSpec> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
struct GetValue<Gaps<TSequence, TSpec> > : Value<Gaps<TSequence, TSpec> >
{};

template <typename TSequence, typename TSpec>
struct GetValue<Gaps<TSequence, TSpec> const> : GetValue<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
struct Position<Gaps<TSequence, TSpec> >
{
    typedef typename Position<TSequence>::Type TSeqPos_;
    typedef typename MakeSigned<TSeqPos_>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Position<Gaps<TSequence, TSpec> const> : Position<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
struct Reference<Gaps<TSequence, TSpec> >
{
    typedef typename Iterator<Gaps<TSequence, TSpec>, Standard>::Type TIterator_;
    typedef Proxy<IteratorProxy<TIterator_> > Type;
};

template <typename TSequence, typename TSpec>
struct Reference<Gaps<TSequence, TSpec> const>
{
    typedef typename Iterator<Gaps<TSequence, TSpec> const, Standard>::Type TIterator_;
    typedef Proxy<IteratorProxy<TIterator_> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
struct Size<Gaps<TSequence, TSpec> >
{
    typedef typename Size<TSequence>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Size<Gaps<TSequence, TSpec> const> : Size<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Source
// ----------------------------------------------------------------------------

// TODO(holtgrew): Switch to Hosted Type interface

template <typename TSequence, typename TSpec>
struct Source<Gaps<TSequence, TSpec> >
{
    typedef TSequence Type;
};

template <typename TSequence, typename TSpec>
struct Source<Gaps<TSequence, TSpec> const> : Source<Gaps<TSequence, TSpec> >
{};

// TODO(holtgrew): Also prefix/suffix/infix? Should work!

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
struct IsSequence<Gaps<TSequence, TSpec> >
{
    typedef True Type;
    static const bool VALUE = true;
};

template <typename TSequence, typename TSpec>
struct IsSequence<Gaps<TSequence, TSpec> const> : IsSequence<Gaps<TSequence, TSpec> >
{};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

// From ContainerConcept, only overwriting documentation here.

/*!
 * @fn Gaps#iter
 * @brief Return an iterator to a specific position in the current clipping.
 *
 * @signature TIterator iter(gaps, viewPos[, tag]);
 *
 * @param[in] gaps    The Gaps object to get an iterator into.
 * @param[in] viewPos View position to get an iterator to (Metafunction: @link ContainerConcept#Position @endlink).
 * @param[in] tag     An optional tag for selecting the iterator type.
 *
 * @return TIterator The resulting iterator.  The type is <tt>Iterator<TGaps, TTag>::Type</tt> where <tt>TTag</tt> is
 *                   the type of <tt>tag</tt>.
 */

// TODO(holtgrew): Adding links to implemented sequence. This should be cleaned up once we have better documentation with concepts.
// ----------------------------------------------------------------------------
// Function setSource()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function createSource()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function clearClipping()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#clearClipping
 * @brief Clear clipping from Gaps objects.
 *
 * @signature void clearClipping(gaps);
 *
 * @param[in,out] gaps Object to clear clipping from.
 */

// ----------------------------------------------------------------------------
// Function clearGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#clearGaps
 * @brief Clear gaps from Gaps objects.
 *
 * @signature void clearGaps(gaps);
 *
 * @param[in,out] gaps Object to clear gaps from.
 */

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// From ContainerConcept, only overwriting documentation here.

/*!
 * @fn Gaps#length
 * @brief Return number of gap and characters between the beginning and the end of the clipping.
 *
 * @signature TSize length(gaps);
 *
 * @param[in] gaps The @link Gaps @endlink object to query for its length.
 *
 * @return TSize The number of gaps and characters between the beginning and the end of the clipping (Metafunction:
 *               @link ContainerConcept#Size @endlink).
 */

// ----------------------------------------------------------------------------
// Function unclippedLength()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#unclippedLength
 * @brief Return length of the gapped sequence without clipping.
 *
 * @signature TSize unclippedLength(gaps);
 *
 * @param[in] gaps The Gaps object to query.
 *
 * @return TSize The result (Metafunction: @link ContainerConcept#Size @endlink).
 */

// ----------------------------------------------------------------------------
// Function toViewPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#toViewPosition
 * @brief Conversion from source (without gaps or clipping) to view position (including gaps and clipping).
 *
 * @signature TPos toViewPosition(gaps, sourcePos);
 *
 * @param[in] gaps      The gaps object to use for translation.
 * @param[in] sourcePos The source position (in the underlying sequence) to translate.
 *
 * @return TPos The resulting position in the view (Metafunction: @link ContainerConcept#Position @endlink).
 */

// ----------------------------------------------------------------------------
// Function toSourcePosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#toSourcePosition
 * @brief Conversion from view position (including gaps and clipping) to source (without gaps or clipping).
 *
 * @signature TPos toSourcePosition(gaps, viewPos);
 *
 * @param[in] gaps      The gaps object to use for translation.
 * @param[in] sourcePos The view position (including gaps and clipping) to translate.
 *
 * @return TPos The resulting position in the underlying sequence (Metafunction: @link ContainerConcept#Position
 *              @endlink).
 */

// ----------------------------------------------------------------------------
// Function isGap()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#isGap
 * @brief Query positions in a Gaps object for being a gap.
 *
 * @signature bool isGap(gaps, viewPos);
 *
 * @param[in] gaps    The Gaps object to query.
 * @param[in] viewPos The view position (including clipping and gaps).
 *
 * @return bool The query result.
 */

template <typename TSequence, typename TSpec, typename TPos>
bool isGap(Gaps<TSequence, TSpec> const & gaps, TPos clippedViewPos)
{
    return isGap(iter(gaps, clippedViewPos, Standard()));
}

// ----------------------------------------------------------------------------
// Function isCharacter()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#isCharacer
 * @brief Query positions in a Gaps object for being a character.
 *
 * @signature bool isGap(gaps, viewPos);
 *
 * @param[in] gaps    The Gaps object to query.
 * @param[in] viewPos The view position (including clipping and gaps).
 *
 * @return bool The query result.
 */

template <typename TSequence, typename TSpec, typename TPos>
bool isCharacter(Gaps<TSequence, TSpec> const & gaps, TPos clippedViewPos)
{
    return isCharacter(iter(gaps, clippedViewPos, Standard()));
}

// ----------------------------------------------------------------------------
// Function insertGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#insertGaps
 * @brief Insert gap characters.
 *
 * @signature void insertGaps(gaps, viewPos, count);
 *
 * @param[in,out] gaps    The Gaps object to insert gaps into.
 * @param[in]     viewPos The view position (including clipping and gaps) to insert gaps at.
 * @param[in]     count   The number of gaps to insert.
 */

// ----------------------------------------------------------------------------
// Function insertGap()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#insertGap
 * @brief Insert one gap character.
 *
 * @signature void insertGap(gaps, viewPos);
 *
 * @param[in,out] gaps    The Gaps object to insert gap into.
 * @param[in]     viewPos The view position (including clipping and gaps) to insert the gap at.
 */

// Forward to removeGaps() which has to be implemented in each subclass.

template <typename TSequence, typename TSpec, typename TPosition>
inline void
insertGap(Gaps<TSequence, TSpec> & gaps, TPosition clippedViewPos)
{
    insertGaps(gaps, clippedViewPos, 1u);
}

// ----------------------------------------------------------------------------
// Function removeGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#removeGaps
 * @brief Remove gaps from a Gaps object.
 *
 * @signature TSize removeGaps(gaps, viewPos, count);
 *
 * @param[in,out] gaps    The gaps object to remove gap characters from.
 * @param[in]     viewPos The view positions to remove gap characters from.
 * @param[in]     count   The number of gap characters to remove.
 *
 * @return TSize The number of gap characters removed (Metafunction: @link ContainerConcept#Size @endlink).
 */

// ----------------------------------------------------------------------------
// Function removeGap()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#removeGap
 * @brief Remove one gap from a Gaps object.
 *
 * @signature TSize removeGap(gaps, viewPos);
 *
 * @param[in,out] gaps    The gaps object to remove one gap character from.
 * @param[in]     viewPos The view positions to remove one gap character from.
 *
 * @return TSize The number of gap characters removed (Metafunction: @link ContainerConcept#Size @endlink).
 */

// Forward to removeGaps() which has to be implemented in each subclass.

template <typename TSequence, typename TSpec, typename TPosition>
inline typename Size<Gaps<TSequence, TSpec> >::Type
removeGap(Gaps<TSequence, TSpec> & gaps, TPosition clippedViewPos)
{
    return removeGaps(gaps, clippedViewPos, 1u);
}

// ----------------------------------------------------------------------------
// Function countGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#countGaps
 * @brief The number of gaps following a position.
 *
 * @signature TSize countGaps(gaps, viewPos);
 *
 * @param[in] gaps    The Gaps object to query.
 * @param[in] viewPos View position (including clipping and gaps) to query at.
 *
 * @return TSize The number of gap characters at <tt>viewPos</tt>  (Metafunction: @link ContainerConcept#Size
 *               @endlink).
 */

template <typename TSequence, typename TSpec, typename TPos>
typename Size<Gaps<TSequence, TSpec> >::Type
countGaps(Gaps<TSequence, TSpec> const & gaps, TPos clippedViewPos)
{
    return countGaps(iter(gaps, clippedViewPos, Standard()));
}

// ----------------------------------------------------------------------------
// Function countLeadingGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#countLeadingGaps
 * @brief The number of leading gaps.
 *
 * @signature TSize countLeadingGaps(gaps);
 *
 * @param[in] gaps    The Gaps object to query.
 *
 * @return TSize The number of leading gap characters  (Metafunction: @link ContainerConcept#Size @endlink).
 */

template <typename TGaps>
inline typename Size<TGaps>::Type
countLeadingGaps(TGaps const & gaps)
{
    return toViewPosition(gaps, 0);
}

// ----------------------------------------------------------------------------
// Function countTrailingGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#countTrailingGaps
 * @brief The number of trailing gaps.
 *
 * @signature TSize countTrailingGaps(gaps);
 *
 * @param[in] gaps    The Gaps object to query.
 *
 * @return TSize The number of trailing gap characters  (Metafunction: @link ContainerConcept#Size @endlink).
 */

template <typename TGaps>
inline typename Size<TGaps>::Type
countTrailingGaps(TGaps const & gaps)
{
    return length(gaps) - toViewPosition(gaps, length(source(gaps)) - 1) - 1;
}

// ----------------------------------------------------------------------------
// Function countCharacters()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#countCharacters
 * @brief The number of characters following a position.
 *
 * @signature TSize countCharacters(gaps, viewPos);
 *
 * @param[in] gaps    The Gaps object to query.
 * @param[in] viewPos View position (including clipping and gaps) to query at.
 *
 * @return TSize The number of non-gaps characters characters at <tt>viewPos</tt> (Metafunction: @link
 *               ContainerConcept#Size @endlink).
 */

template <typename TSequence, typename TSpec, typename TPos>
typename Size<Gaps<TSequence, TSpec> >::Type
countCharacters(Gaps<TSequence, TSpec> const & gaps, TPos clippedViewPos)
{
    return countCharacters(iter(gaps, clippedViewPos, Standard()));
}

// ----------------------------------------------------------------------------
// Function setClippedBeginPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#setClippedBeginPosition
 * @brief Set the begin position of the clipping.
 *
 * @signature void setClippedBeginPosition(gaps, unclippedViewPos);
 *
 * @param[in,out] gaps             The Gaps object to set the clipping begin position of.
 * @param[in]     unclippedViewPos View position (including gaps but excluding clipping) to set the clipping begin to.
 */

// ----------------------------------------------------------------------------
// Function setClippedEndPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#setClippedEndPosition
 * @brief Set the end position of the clipping.
 *
 * @signature void setClippedEndPosition(gaps, unclippedViewPos);
 *
 * @param[in,out] gaps             The Gaps object to set the clipping end position of.
 * @param[in]     unclippedViewPos View position (including gaps but excluding clipping) to set the clipping end to.
 */

// ----------------------------------------------------------------------------
// Function clippedBeginPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#clippedBeginPosition
 * @brief Return begin position of the clipping.
 *
 * @signature TPos clippedBeginPosition(gaps);
 *
 * @param[in] gaps             The Gaps object to query.
 *
 * @return TPos The begin position of the unclipped view  (Metafunction: @link ContainerConcept#Position @endlink).
 *
 * @section Example
 *
 * In the following gaps configuration, the result of <tt>clippedBeginPosition(gaps)</tt> is 1.
 *
 * @code{.txt}
 * clipping                   [     )
 *   (half-open interval)
 *
 * gapped sequence:          X--XXX-XX-
 *
 * source position:          0111234456
 * unclipped view position:  0123456789
 * clipped view position:     0123456
 * @endcode
 */

// ----------------------------------------------------------------------------
// Function clippedEndPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#clippedEndPosition
 * @brief Return end position of the clipping.
 *
 * @signature TPos clippedEndPosition(gaps);
 *
 * @param[in] gaps             The Gaps object to query.
 *
 * @return TPos The end position of the unclipped view  (Metafunction: @link ContainerConcept#Position @endlink).
 *
 * @section Example
 *
 * In the following gaps configuration, the result of <tt>clippedEndPosition(gaps)</tt> is 7.
 *
 * @code{.txt}
 * clipping                   [     )
 *   (half-open interval)
 *
 * gapped sequence:          X--XXX-XX-
 *
 * source position:          0111234456
 * unclipped view position:  0123456789
 * clipped view position:     0123456
 * @endcode
 */

// ----------------------------------------------------------------------------
// Function setBeginPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#setBeginPosition
 * @brief Set the begin position of the clipped gapped sequence, given a source position.
 *
 * @signature void setBeginPosition(gaps, sourcePos);
 *
 * @param[in,out] gaps      The Gaps object to set the begin position in.
 * @param[in]     sourcePos Position in the underlying sequence to set clipping to.
 */

// ----------------------------------------------------------------------------
// Function setEndPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#setEndPosition
 * @brief Set the end position of the clipped gapped sequence, given a source position.
 *
 * @signature void setEndPosition(gaps, sourcePos);
 *
 * @param[in,out] gaps      The Gaps object to set the end position in.
 * @param[in]     sourcePos Position in the underlying sequence to set clipping to.
 */

// ----------------------------------------------------------------------------
// Function beginPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#beginPosition
 * @brief Return the clipping begin position as a source position.
 *
 * @signature TPos beginPosition(gaps);
 *
 * @param[in] gaps The Gaps object to query.
 *
 * @return TPos The clipping begin position in the source (Metafunction: @link ContainerConcept#Position @endlink).
 *
 * @section Example
 *
 * In the following gaps configuration, the result of <tt>beginPosition(gaps)</tt> is $1$.  The clipping starts in a
 * gap and the source position of the first non-gap character right of the clipping begin has source position 1.
 *
 * @code{.txt}
 * clipping                   [     )
 *   (half-open interval)
 *
 * gapped sequence:          X--XXX-XX-
 *
 * source position:          0111234456
 * unclipped view position:  0123456789
 * clipped view position:     0123456
 * @endcode
 */

// ----------------------------------------------------------------------------
// Function endPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#endPosition
 * @brief Return the clipping end position as a source position.
 *
 * @signature TPos endPosition(gaps);
 *
 * @param[in] gaps The Gaps object to query for the end position as a source position.
 *
 * @return TPos The end position as a source position (Metafunction: @link ContainerConcept#Position @endlink).
 *
 * @section Example
 *
 * In the following gaps configuration, the result of <tt>endPositioN(gaps)</tt> is 4.
 *
 * @code{.txt}
 * clipping                   [     )
 *   (half-open interval)
 *
 * gapped sequence:          X--XXX-XX-
 *
 * source position:          0111234456
 * unclipped view position:  0123456789
 * clipped view position:     0123456
 * @endcode
 */

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSource, typename TSpec>
inline void
write(TTarget & target,
      Gaps<TSource, TSpec> const & source)
{
    // Print gaps row
    typedef typename Iterator<Gaps<TSource, TSpec> const>::Type TIter;
    TIter begin_ = begin(source);
    TIter end_ = end(source);
    for (; begin_ != end_; ++begin_) {
        if (isGap(begin_))
            writeValue(target, gapValue<char>());
        else
            writeValue(target, convert<char>(getValue(begin_)));
    }
}

// ----------------------------------------------------------------------------
// Function operator<<()                                      [stream operator]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document appropriately.

template <typename TTarget, typename TSource, typename TSpec>
inline TTarget &
operator<<(TTarget & target, Gaps<TSource, TSpec> const & gaps)
{
    typename DirectionIterator<TTarget, Output>::Type it = directionIterator(target, Output());
    write(it, gaps);
    return target;
}

// ----------------------------------------------------------------------------
// Function _pumpTraceToGaps()
// ----------------------------------------------------------------------------

// Internal function for converting AlignTrace<> objects into alignments in two Gaps objects.  Note that the traceback
// in the trace is stored in reverse, from back to front.  We insert the gaps in descending order of their position.
// The reason is that Gaps<> objects store the gaps in ascending order of coordinates in String<> objects and inserting
// at the end is in O(1) while inserting in the front is O(n).

template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TSize>
void _pumpTraceToGaps(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                      Gaps<TSequenceV, TGapsSpecV> & gapsV,
                      AlignTraceback<TSize> const & trace)
{
    typedef Gaps<TSequenceH, TGapsSpecH> TGapsH;
    typedef typename Iterator<TGapsH, Standard>::Type TGapsHIter;

    typedef Gaps<TSequenceV, TGapsSpecV> TGapsV;
    typedef typename Iterator<TGapsV, Standard>::Type TGapsVIter;

    // TODO(holtgrew): I don't understand the following.  Originally, this function used Align objects, but I did not understand it there either.
    // TODO(rausch): Pump trace into align_ (note: this is relatively slow code here. it could be improved if specialized to the Align Specs).
    clearGaps(gapsH);
    clearClipping(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsV);

    TSize i = length(trace.sizes);  // Scan trace backwards.
    TGapsHIter itH = begin(gapsH);
    TGapsVIter itV = begin(gapsV);
    while (i > 0)
    {
        --i;
        TSize size = trace.sizes[i];
        switch ((int) trace.tvs[i])
        {
        case 1:  // Go horizontal.
            insertGaps(itV, size);
            break;

        case 2:  // Go vertical.
            insertGaps(itH, size);
            break;
        }
        goFurther(itH, size);
        goFurther(itV, size);
    }
}

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document TSource via metafunctio.

/*!
 * @fn Gaps#source
 * @brief Return underlying object.
 *
 * @signature TSource source(gaps);
 *
 * @param[in] gaps The Gaps object to return the underling sequence for.
 *
 * @return TSource Reference to the source of the Gaps.
 */

// ----------------------------------------------------------------------------
// Function sourceSegment()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Rename/remove?

// We need some forwards for function sourceSegment().

template <typename TSequence, typename TSpec>
inline typename Position<Gaps<TSequence, TSpec> >::Type clippedBeginPosition(Gaps<TSequence, TSpec> const & gaps);
template <typename TSequence, typename TSpec>
inline typename Position<Gaps<TSequence, TSpec> >::Type clippedEndPosition(Gaps<TSequence, TSpec> const & gaps);
template <typename TSequence, typename TSpec, typename TPosition>
inline typename Position<TSequence>::Type
toSourcePosition(Gaps<TSequence, TSpec> const & gaps, TPosition clippedViewPos);

template <typename TSequence, typename TSpec>
inline typename Infix<TSequence>::Type
sourceSegment(Gaps<TSequence, TSpec> const & gaps)
{
    return infix(source(gaps), toSourcePosition(gaps, clippedBeginPosition(gaps)), toSourcePosition(gaps, clippedEndPosition(gaps)));
}

template <typename TSequence, typename TSpec>
inline typename Infix<TSequence>::Type
sourceSegment(Gaps<TSequence, TSpec> & gaps)
{
    return infix(source(gaps), toSourcePosition(gaps, clippedBeginPosition(gaps)), toSourcePosition(gaps, clippedEndPosition(gaps)));
}

// ----------------------------------------------------------------------------
// Function assignSource()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#assignSource
 * @brief Assign the source of a gaps object, copying data.
 *
 * @signature void assignSource(gaps, seq);
 *
 * @param[in,out] gaps The Gaps object to assign the source of.
 * @param[in]     seq  The @link ContainerConcept sequence @endlink to assign to the underlying string.
 */

// TOOD(holtgrew): Switch to Hosted Type?

template <typename TSequence, typename TSpec, typename TValue>
inline void
assignSource(Gaps<TSequence, TSpec> & gaps, TValue const & value)
{
    assign(source(gaps), value);
}

// ----------------------------------------------------------------------------
// Function setSource()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function copyGaps()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#copyGaps
 * @brief Copy gaps from one Gaps object to another (in the clipped view of both argumetns).
 *
 * The user is resposible for ensuring that the gaps are over sequences of same length and appropriate clipping.
 *
 * @signature void copyGaps(dest, source);
 *
 * @param[in,out] dest   The destination Gaps object (appropriate clipping, no gaps).
 * @param[in]     source The source Gaps object.
 */

template <typename TDestSource, typename TDestSpec, typename TSourceSource, typename TSourceSpec>
void copyGaps(Gaps<TDestSource, TDestSpec> & dest, Gaps<TSourceSource, TSourceSpec> const & source)
{
    typedef Gaps<TDestSource, TDestSpec> TLhs;
    typedef typename Iterator<TLhs, Standard>::Type TLhsIter;
    typedef Gaps<TSourceSource, TSourceSpec> const TRhs;
    typedef typename Iterator<TRhs, Standard>::Type TRhsIter;

    TLhsIter lhsIt = begin(dest, Standard());
    //TLhsIter lhsItEnd = end(dest, Standard());
    TRhsIter rhsIt = begin(source, Standard());
    TRhsIter rhsItEnd = end(source, Standard());

    for (unsigned num = 0; rhsIt != rhsItEnd; lhsIt += num, rhsIt += num)
    {
        if (isGap(rhsIt))
        {
            num = countGaps(rhsIt);
            insertGaps(lhsIt, num);
        }
        else
        {
            num = countCharacters(rhsIt);
        }

        SEQAN_ASSERT_NOT((lhsIt == end(dest, Standard())) && num > 0);
    }
}

// ----------------------------------------------------------------------------
// Function copyClipping()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#copyClipping
 * @brief Copy clipping information from one Gaps object to another.
 *
 * @signature void copyClipping(dest, source);
 *
 * @param[in,out] dest   The destination Gaps object.
 * @param[in]     source The source Gaps object.
 */

template <typename TDestSource, typename TDestSpec, typename TSourceSource, typename TSourceSpec>
void copyClipping(Gaps<TDestSource, TDestSpec> & dest, Gaps<TSourceSource, TSourceSpec> const & source)
{
    setClippedBeginPosition(dest, clippedBeginPosition(source));
    setClippedEndPosition(dest, clippedEndPosition(source));
}

// ----------------------------------------------------------------------------
// Function clipSemiGlobal()
// ----------------------------------------------------------------------------

/*!
 * @fn Gaps#clipSemiGlobal
 * @brief Clip the Gaps objects to reflect a semi-global alignment.
 *
 * Leading and trailing gaps are clipped in the local Gaps object. The global Gaps object is updated accordingly.
 *
 * @signature void clipSemiGlobal(global, local);
 *
 * @param[in,out] global The global Gaps object.
 * @param[in,out] local  The local Gaps object.
 */

template <typename TGlobalGaps, typename TLocalGaps>
inline void clipSemiGlobal(TGlobalGaps & global, TLocalGaps & local)
{
    typedef typename Size<TLocalGaps>::Type  TGapsSize;

    TGapsSize leadingGaps = countLeadingGaps(local);
    TGapsSize trailingGaps = countTrailingGaps(local);
    TGapsSize globalLenght = length(global);
    TGapsSize localLength = length(local);

    setClippedBeginPosition(global, leadingGaps);
    setClippedBeginPosition(local, leadingGaps);
    setClippedEndPosition(global, globalLenght - trailingGaps);
    setClippedEndPosition(local, localLength - trailingGaps);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline void clearGaps(Gaps<TSequence, TSpec> & gaps);
template <typename TSequence, typename TSpec>
inline void clearClipping(Gaps<TSequence, TSpec> & gaps);

template <typename TSequence, typename TSpec>
inline void clear(Gaps<TSequence, TSpec> & gaps)
{
    clearGaps(gaps);
    clearClipping(gaps);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline bool operator==(Gaps<TSequence, TSpec> const & lhs,
                       Gaps<TSequence, TSpec> const & rhs)
{
    typename Comparator<Gaps<TSequence, TSpec> >::Type lex(lhs, rhs);
    return isEqual(lex);
}

template <typename TSequence, typename TSpec, typename TRightHandSide>
inline bool operator==(Gaps<TSequence, TSpec> const & lhs,
                       TRightHandSide const & rhs)
{
    typename Comparator<Gaps<TSequence, TSpec> >::Type lex(lhs, rhs);
    return isEqual(lex);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline bool operator!=(Gaps<TSequence, TSpec> const & lhs,
                       Gaps<TSequence, TSpec> const & rhs)
{
    return !(lhs == rhs);
}


}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_
