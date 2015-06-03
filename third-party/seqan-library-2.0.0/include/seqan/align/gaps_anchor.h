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
// Implementation of Gaps class using pairs of gapped sequence and source
// sequence position as anchors.
// ==========================================================================

// TODO(holtgrew): Clipping in leading and trailing gaps is not possible right now. Dave and I have to discuss this further.
// TODO(holtgrew): Also, inserting gaps in the front changes the clipped begin position which is unexpected.

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ANCHOR_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ANCHOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

/*!
 * @tag GapsSpecTag#AnchorGaps
 * @headerfile <seqan/align.h>
 * @brief Tag for the Anchor Gaps specialization.
 *
 * @signature template <typename TGapAnchors>
 *            struct AnchorGaps;
 */

template <typename TGapAnchors>
struct AnchorGaps;

template <typename TSequence, typename TGapAnchors>
inline void _reinitAnchorGaps(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Specialization AnchorGaps
// ----------------------------------------------------------------------------

/*!
 * @class AnchorGaps
 * @extends Gaps
 * @headerfile <seqan/align.h>
 * @brief Stores gaps as anchors of the first characters behind gaps.
 *
 * @signature template <typename TSource, typename TGapAnchors = String<GapAnchor<unsigned> > >
 *            class Gaps<TSource, AnchorGaps<TGapAnchors> >;
 *
 * @tparam TSource     The type of the underling sequence.
 * @tparam TGapAnchors The type of the string of @link GapAnchor @endlink objects.
 */

/*!
 * @fn AnchorGaps::Gaps
 * @brief Constructor
 *
 * @signature Gaps::Gaps([other]);
 * @signature Gaps::Gaps(source[, anchors]);
 * @signature Gaps::Gaps(anchors);
 *
 * @param[in] other   Another @link AnchorGaps @endlink object to copy from.
 * @param[in] source  The underling sequence to construct the Gaps object from.
 * @param[in] anchors The string of anchors to construct with.
 *
 * An AnchorGaps object has a default constructor, can be constructed from the underlying source, and/or a string of
 * gap anchors.
 */


template <typename TGapAnchors = String<GapAnchor<unsigned> > >
struct AnchorGaps
{};

template <typename TSource, typename TGapAnchors>
class Gaps<TSource, AnchorGaps<TGapAnchors> >
{
public:
    // -----------------------------------------------------------------------
    // Internal Typedefs
    // -----------------------------------------------------------------------

    typedef typename Value<TGapAnchors>::Type TGapAnchor_;
    typedef typename Position<TGapAnchor_>::Type TViewPosition_;
    typedef typename Position<Gaps>::Type      TPosition_;
    typedef typename Value<Gaps>::Type         TValue_;

    // -----------------------------------------------------------------------
    // Member Variables
    // -----------------------------------------------------------------------

    Holder<TSource>     data_source;
    Holder<TGapAnchors> data_gaps;
    int                 data_cutBegin;      // number of gap positions cut from the beginning
    int                 data_cutEnd;        // number of gap positions cut from the end
    int                 data_viewCutBegin;  // how many alignment chars should be clipped at the beginning (can be negative too)
    int                 data_viewCutEnd;    // how ...                                           end ...

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    Gaps() :
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    Gaps(TSource & source) :
        data_source(source),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    Gaps(TGapAnchors & anchors) :
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    Gaps(TGapAnchors const & anchors) :
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    // Note: We need the variants with the first parameter "TSource const &" here because TSource can be a Segment which
    // is often given as a temporary.

    Gaps(TSource & source, TGapAnchors & anchors) :
        data_source(source),
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    Gaps(TSource & source, TGapAnchors const & anchors) :
        data_source(source),
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    // TODO(holtgrew): These constructors are only here because of const-Holder issues.

    template <typename TSource2>
    Gaps(TSource2 & source, TGapAnchors & anchors) :
        data_source(source),
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    template <typename TSource2>
    Gaps(TSource2 & source, TGapAnchors const & anchors) :
        data_source(source),
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    template <typename TSource2>
    Gaps(TSource2 const & source, TGapAnchors & anchors) :
        data_source(source),
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    template <typename TSource2>
    Gaps(TSource2 const & source, TGapAnchors const & anchors) :
        data_source(source),
        data_gaps(anchors),
        data_cutBegin(0),
        data_cutEnd(0),
        data_viewCutBegin(0),
        data_viewCutEnd(0)
    {
    }

    // -----------------------------------------------------------------------
    // Array Subscript Operator
    // -----------------------------------------------------------------------

    inline TValue_
    operator[](TPosition_ clippedViewPos) const
    {
        return value(*this, clippedViewPos);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function swap()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors>
void swap(Gaps<TSequence, AnchorGaps<TGapAnchors> > & lhs,
          Gaps<TSequence, AnchorGaps<TGapAnchors> > & rhs)
{
    swap(lhs.data_source, rhs.data_source);
    swap(lhs.data_gaps, rhs.data_gaps);

    std::swap(lhs.data_cutBegin, rhs.data_cutBegin);
    std::swap(lhs.data_cutEnd, rhs.data_cutEnd);
    std::swap(lhs.data_viewCutBegin, rhs.data_viewCutBegin);
    std::swap(lhs.data_viewCutEnd, rhs.data_viewCutEnd);
}

// ----------------------------------------------------------------------------
// Helper Function _reinitAnchorGaps()
// ----------------------------------------------------------------------------

// Reset the array gaps DS such that represents the ungapped sequence.

template <typename TSequence, typename TGapAnchors>
inline void _reinitAnchorGaps(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps)
{
    clear(value(gaps.data_gaps));
    gaps.data_cutBegin = 0;
    gaps.data_cutEnd = 0;
    gaps.data_viewCutBegin = 0;
    gaps.data_viewCutEnd = 0;
}

// ----------------------------------------------------------------------------
// Function _dataSource()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors>
inline Holder<TSource> &
_dataSource(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
    return me.data_source;
}

template <typename TSource, typename TGapAnchors>
inline Holder<TSource> const &
_dataSource(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
    return me.data_source;
}

// ----------------------------------------------------------------------------
// Function _assignSourceLength()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSource, typename TGapAnchors>
inline void
_assignSourceLength(TSize & size, Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
    if (IsSameType<TSource, Nothing>::VALUE)
        size = maxValue<TSize>() / 2;
    else
        size = length(value(me.data_source));
}

// ----------------------------------------------------------------------------
// Function _dataAnchors()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors>
inline TGapAnchors &
_dataAnchors(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
    return value(me.data_gaps);
}

template <typename TSource, typename TGapAnchors>
inline TGapAnchors const &
_dataAnchors(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
    return value(const_cast<Holder<TGapAnchors> &>(me.data_gaps));
}

// ----------------------------------------------------------------------------
// Function _getAnchor()
// ----------------------------------------------------------------------------

template <typename TAnchor, typename TSource, typename TGapAnchors, typename TIdx>
inline void
_getAnchor(TAnchor & anchor, Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, TIdx idx)
{
    if (idx > (TIdx)length(_dataAnchors(me)))
    {
        _assignSourceLength(anchor.seqPos, me);
        if (empty(_dataAnchors(me)) && idx == 1)
            anchor.gapPos = anchor.seqPos;
        else
        {
            // for the sick case that an anchor seq position is beyond the sequence end
            if (!empty(_dataAnchors(me)))
            {
                // if there is no sequence but anchors -> assume infinite sequence
                if (anchor.seqPos == 0)
                    anchor.seqPos = maxValue(anchor.gapPos);
                // if the sequence has a length > 0, but there is an anchor behind the end
                // -> elongate sequence
                else if ((__int64)anchor.seqPos < (__int64)back(_dataAnchors(me)).seqPos)
                    anchor.seqPos = back(_dataAnchors(me)).seqPos;
            }
            anchor.gapPos = maxValue(anchor.gapPos);
        }
    }
    else if (idx > 0)
        anchor = _dataAnchors(me)[idx - 1];
    else
    {
        anchor.seqPos = 0;
        if (idx == 0)
            anchor.gapPos = 0;
        else
            anchor.gapPos = minValue(anchor.gapPos);
    }
}

// ----------------------------------------------------------------------------
// Function _unclippedLength()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors>
inline typename Size<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
_unclippedLength(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
    typedef typename Value<TGapAnchors>::Type TAnchor;
    typedef typename Size<TAnchor>::Type TSize;
    TSize len;
    _assignSourceLength(len, me);
    if (!empty(_dataAnchors(me)))
    {
        TAnchor const & last = back(_dataAnchors(me));
        len += last.gapPos - last.seqPos;
    }
    return len - (me.data_cutBegin + me.data_cutEnd);
}

// ----------------------------------------------------------------------------
// Function unclippedLength()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors>
inline typename Size<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
unclippedLength(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
    // TODO(holtgrew): Merge into one public function.
    return _unclippedLength(me);
}

// ----------------------------------------------------------------------------
// Function clearClipping()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors>
inline void
clearClipping(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps)
{
    gaps.data_viewCutBegin = 0;
    gaps.data_viewCutEnd = 0;
}

// ----------------------------------------------------------------------------
// Function clearGaps()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors>
inline void
clearGaps(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps)
{
    _reinitAnchorGaps(gaps);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors>
inline void
clear(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps)
{
    clear(gaps.data_source); // clear source holder

    // clear gaps, but on holder level
    clear(gaps.data_gaps);
    gaps.data_cutBegin = 0;
    gaps.data_cutEnd = 0;
    gaps.data_viewCutBegin = 0;
    gaps.data_viewCutEnd = 0;
}

// ----------------------------------------------------------------------------
// Function isGap()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition>
inline bool
isGap(Gaps<TSequence, AnchorGaps<TGapAnchors> > const & gaps, TPosition clippedViewPos)
{
    // TODO(holtgrew): Implement without iterator?
    return isGap(iter(gaps, clippedViewPos));
}

// ----------------------------------------------------------------------------
// Function insertGaps()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition, typename TCount>
inline void
insertGaps(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps, TPosition clippedViewPos, TCount count)
{
    // TODO(holtgrew): Implement without iterator?
    typedef Gaps<TSequence, AnchorGaps<TGapAnchors> > TGaps;
    typedef typename Iterator<TGaps>::Type TIter;

    TIter it = iter(gaps, clippedViewPos);
    insertGaps(it, count);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition>
inline typename Value<Gaps<TSequence, AnchorGaps<TGapAnchors> > >::Type
value(Gaps<TSequence, AnchorGaps<TGapAnchors> > const & gaps, TPosition clippedViewPos)
{
    // TODO(holtgrew): Implement without iterator?
    typedef Gaps<TSequence, AnchorGaps<TGapAnchors> > TGaps;
    typedef typename Iterator<TGaps const>::Type TIter;

    TIter it = iter(gaps, clippedViewPos);
    if (isGap(it))
        return '-';
    return *it;
}

// ----------------------------------------------------------------------------
// Function removeGaps()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition, typename TCount>
inline typename Size<Gaps<TSequence, AnchorGaps<TGapAnchors> > >::Type
removeGaps(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps, TPosition clippedViewPos, TCount count)
{
    // TODO(holtgrew): Implement without iterator?
    return removeGaps(iter(gaps, clippedViewPos), count);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We would rather like to have the const version only :(
template <typename TSource, typename TGapAnchors>
inline typename Size<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
length(Gaps<TSource, AnchorGaps<TGapAnchors> > & me)
{
    return _unclippedLength(me) - (me.data_viewCutBegin + me.data_viewCutEnd);
}

template <typename TSource, typename TGapAnchors>
inline typename Size<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
length(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me)
{
    return _unclippedLength(me) - (me.data_viewCutBegin + me.data_viewCutEnd);
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, Standard)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type(me);
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, Standard)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type(me);
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, Rooted)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type(me);
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type
begin(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, Rooted)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type(me);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
end(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, Standard)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type(me, length(me));
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type
end(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, Standard)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type(me, length(me));
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type
end(Gaps<TSource, AnchorGaps<TGapAnchors> > & me, Rooted)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > >::Type(me, length(me));
}

template <typename TSource, typename TGapAnchors>
inline typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type
end(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, Rooted)
{
    SEQAN_CHECKPOINT
    return typename Iterator<Gaps<TSource, AnchorGaps<TGapAnchors> > const>::Type(me, length(me));
}

// ----------------------------------------------------------------------------
// Function beginPosition()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We would rather like to have the const version only :(
template <typename TSource, typename TGapAnchors>
inline typename Position<TSource>::Type
beginPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & gaps)
{
    return toSourcePosition(gaps, 0);
}

template <typename TSource, typename TGapAnchors>
inline typename Position<TSource>::Type
beginPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > const & gaps)
{
    return toSourcePosition(gaps, 0);
}

// ----------------------------------------------------------------------------
// Function endPosition()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We would rather like to have the const version only :(
template <typename TSource, typename TGapAnchors>
inline typename Position<TSource>::Type
endPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & gaps)
{
    return toSourcePosition(gaps, _unclippedLength(gaps) - (gaps.data_viewCutEnd + gaps.data_viewCutBegin));
}

template <typename TSource, typename TGapAnchors>
inline typename Position<TSource>::Type
endPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > const & gaps)
{
    return toSourcePosition(gaps, _unclippedLength(gaps) - (gaps.data_viewCutEnd + gaps.data_viewCutBegin));
}

// ----------------------------------------------------------------------------
// Function setBeginPosition()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors, typename TPosition>
inline void
setBeginPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & gaps, TPosition sourcePosition)
{
    setClippedBeginPosition(gaps, toViewPosition(gaps, sourcePosition) + clippedBeginPosition(gaps));
}

// ----------------------------------------------------------------------------
// Function setEndPosition()
// ----------------------------------------------------------------------------

template <typename TSource, typename TGapAnchors, typename TPosition>
inline void
setEndPosition(Gaps<TSource, AnchorGaps<TGapAnchors> > & gaps, TPosition sourcePosition)
{
    setClippedEndPosition(gaps, toViewPosition(gaps, sourcePosition) + clippedBeginPosition(gaps));
}

// ----------------------------------------------------------------------------
// Function createSource()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove? Switch to Hosted Type Interface?

template <typename TSequence, typename TGapAnchors>
inline void createSource(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps)
{
    create(gaps.data_source);
}

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchor>
inline typename Source<Gaps<TSequence, AnchorGaps<TGapAnchor> > >::Type &
source(Gaps<TSequence, AnchorGaps<TGapAnchor> > const & gaps)
{
    return value(gaps.data_source);
}

// ----------------------------------------------------------------------------
// Function setSource()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchor>
inline void
setSource(Gaps<TSequence, AnchorGaps<TGapAnchor> > & gaps, TSequence & source)
{
    setValue(gaps.data_source, source);
}

template <typename TSequence, typename TGapAnchor>
inline void
setSource(Gaps<TSequence const, AnchorGaps<TGapAnchor> > & gaps, TSequence & source)
{
    setValue(gaps.data_source, source);
}

// ----------------------------------------------------------------------------
// Function assignSource()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchor, typename TSequence2>
inline void
assignSource(Gaps<TSequence, AnchorGaps<TGapAnchor> > & gaps, TSequence2 const & source)
{
    value(gaps.data_source) = source;
}

// ----------------------------------------------------------------------------
// Function positionGapToSeq()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Deprecate in favour of toViewPosition/toSourcePosition?

/*!
 * @fn AnchorGaps#positionGapsToSeq
 * @brief Convert from gap space in the global alignment to the sequence space on the reference.
 *
 * @signature TPos positionGapToSeq(gaps, pos);
 *
 * @param[in] gaps Contig AnchorGaps (e.g. from FragmentStore).
 * @param[in] pos  Position in gap space.
 *
 * @return TPos Position in sequence space (Metafunction: @link ContainerConcept#Position @endlink).
 *
 * See the example below to construct the Gaps ojbect.  Note that this construction is fast since it ionly a thing wrapper
 * around underlying objects.
 *
 * @section Examples
 *
 * Convert from gap space to positions pace when the contig required to be loaded. * Converts position aligned read with
 * index <tt>idx</tt> in the aligned read store.
 *
 * @code{.cpp}
 * typedef typename TFragmentStore::TContigStore                        TContigStore;
 * typedef typename Value<TContigStore>::Type                           TContig;
 * typedef typename TFragmentStore::TContigSeq                          TContigSeq;
 * typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
 *
 * typedef typename TFragmentStore::TAlignedReadStore                   TAlignedReadStore;
 * typedef typename Value<TAlignedReadStore>::Type                      TAlignedRead;
 * typedef typename TAlignedRead::TPos                                  TAlignedReadPos;
 *
 * unsigned contigId = alignedReadStore[idx].contigId;
 * TContigGaps contigGaps(contigStore[contigId].seq, contigStore[contigId].gaps);
 * TAlignedRead const & alignedRead = alignedReadStore[idx];
 * // Translate end position from aligned read record to sequence space in reference.
 * TAlignedReadPos endPos = positionGapToSeq(contigGaps, alignedRead.endPos);
 * @endcode
 *
 * Convert from gap space to position space when the contigs are not required.
 * Converts position aligned read with index $idx$ in the aligned read store.
 *
 * @code{.cpp}
 * typedef typename TFragmentStore::TContigStore                        TContigStore;
 * typedef typename Value<TContigStore>::Type                           TContig;
 * typedef Gaps<Nothing, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
 *
 * typedef typename TFragmentStore::TAlignedReadStore                   TAlignedReadStore;
 * typedef typename Value<TAlignedReadStore>::Type                      TAlignedRead;
 * typedef typename TAlignedRead::TPos                                  TAlignedReadPos;
 *
 * unsigned contigId = alignedReadStore[idx].contigId;
 * TContigGaps contigGaps(Nothing(), contigStore[contigId].gaps);
 * TAlignedRead const & alignedRead = alignedReadStore[idx];
 * // Translate end position from aligned read record to sequence space in reference.
 * TAlignedReadPos endPos = positionGapToSeq(contigGaps, alignedRead.endPos);
 * @endcode
 */

template <typename TSource, typename TGapAnchors, typename TPosition>
inline TPosition
positionGapToSeq(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, TPosition pos)
{
    typedef typename Position<typename Value<TGapAnchors>::Type>::Type TAnchorPos;

    GapAnchor<__int64> prevAnchor, nextAnchor;
    TPosition           seqPos;
    int                 anchorIdx;

    if (isNegative(pos))
        anchorIdx = -1;
    else
    {
        TGapAnchors const & anchors = _dataAnchors(me);
        TAnchorPos seqLength;
        _assignSourceLength(seqLength, me);
        if (!empty(anchors))
        {
            anchorIdx = upperBoundGapAnchor(anchors, pos, SortGapPos()) - begin(anchors, Standard());
            if (anchorIdx < (int)length(anchors))
                if (anchors[anchorIdx].gapPos == (TAnchorPos)pos && anchors[anchorIdx].seqPos != seqLength)
                    ++anchorIdx;
        }
        else
            anchorIdx = ((TAnchorPos)pos < seqLength) ? 0 : 1;
    }
    _getAnchor(prevAnchor, me, anchorIdx);
    _getAnchor(nextAnchor, me, anchorIdx + 1);

    if (nextAnchor.seqPos - prevAnchor.seqPos > (int)pos - prevAnchor.gapPos)
        seqPos = prevAnchor.seqPos + (pos - prevAnchor.gapPos);
    else
        seqPos = nextAnchor.seqPos;
    return seqPos;
}

// ----------------------------------------------------------------------------
// Function positionSeqToGap()
// ----------------------------------------------------------------------------

/*!
 * @fn AnchorGaps#positionSeqToGap
 * @brief Convert from sequence space on the reference to gap space in the global alignment.
 *
 * @signature TPos positionSeqToGap(gaps, pos);
 *
 * @param[in] gaps The AnchorGaps object to use for the translation.
 * @param[in] pos  The gap space position to conver to sequence space.
 *
 * @return TPos The resulting position in sequence space (Metafunction: @link ContainerConcept#Position @endlink).
 *
 * See the example below to construct the gaps object.  Note that this construction is fast since it is only a thin
 * wrapper around underlying objects.
 *
 * @section Example
 *
 * Convert from gap space to position space on contig $contigId$ when the contigs required to be loaded.
 *
 * @code{.cpp}
 * typedef typename TFragmentStore::TContigStore                        TContigStore;
 * typedef typename Value<TContigStore>::Type                           TContig;
 * typedef typename TFragmentStore::TContigSeq                          TContigSeq;
 * typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
 *
 * TContigGaps contigGaps(contigStore[contigId].seq, contigStore[contigId].gaps);
 * TAlignedReadPos pos = positionGapToSeq(contigGaps, 33);
 * @endcode
 *
 * Convert from gap space to position space on contig $contigId$ when the contigs are not required.
 *
 * @code{.cpp}
 * typedef typename TFragmentStore::TContigStore                        TContigStore;
 * typedef typename Value<TContigStore>::Type                           TContig;
 * typedef Gaps<Nothing, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
 *
 * TContigGaps contigGaps(Nothing(), contigStore[contigId].gaps);
 * TAlignedReadPos endPos = positionGapToSeq(contigGaps, 33);
 * @endcode
 */

template <typename TSource, typename TGapAnchors, typename TPosition>
inline TPosition
positionSeqToGap(Gaps<TSource, AnchorGaps<TGapAnchors> > const & me, TPosition pos)
{
    typedef typename Position<typename Value<TGapAnchors>::Type>::Type TAnchorPos;

    GapAnchor<__int64>  prevAnchor, nextAnchor;
    TPosition           gapPos;
    int                 anchorIdx;

    if (isNegative(pos))
        anchorIdx = -1;
    else
    {
        TGapAnchors const & anchors = _dataAnchors(me);
        TAnchorPos seqLength;
        _assignSourceLength(seqLength, me);
        if (!empty(anchors))
        {
            anchorIdx = upperBoundGapAnchor(anchors, pos, SortSeqPos()) - begin(anchors, Standard());
            if (anchorIdx < (int)length(anchors))
                if (anchors[anchorIdx].seqPos == (TAnchorPos)pos)
                    ++anchorIdx;
        }
        else
            anchorIdx = ((TAnchorPos)pos < seqLength) ? 0 : 1;
    }
    _getAnchor(prevAnchor, me, anchorIdx);
    _getAnchor(nextAnchor, me, anchorIdx + 1);

    if (nextAnchor.gapPos - prevAnchor.gapPos > (int)pos - prevAnchor.seqPos)
        gapPos = prevAnchor.gapPos + (pos - prevAnchor.seqPos);
    else
        gapPos = nextAnchor.gapPos;
    return gapPos;
}

// ----------------------------------------------------------------------------
// Function toViewPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition>
inline typename Position<Gaps<TSequence, AnchorGaps<TGapAnchors> > >::Type
toViewPosition(Gaps<TSequence, AnchorGaps<TGapAnchors> > const & gaps, TPosition sourcePosition)
{
    return positionSeqToGap(gaps, sourcePosition) - gaps.data_viewCutBegin - gaps.data_cutBegin;
}

// ----------------------------------------------------------------------------
// Function toSourcePosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition>
inline typename Position<TSequence>::Type
toSourcePosition(Gaps<TSequence, AnchorGaps<TGapAnchors> > const & gaps, TPosition clippedViewPos)
{
    // TODO(weese): possibly change positionGapToSeq interface to consider a different zero
    // shifted by data_cutBegin
    return positionGapToSeq(gaps, clippedViewPos + gaps.data_viewCutBegin + gaps.data_cutBegin);
}

// ----------------------------------------------------------------------------
// Function setClippedBeginPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition>
inline void
setClippedBeginPosition(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps, TPosition unclippedViewPosition)
{
    gaps.data_viewCutBegin = unclippedViewPosition;
}

// ----------------------------------------------------------------------------
// Function setClippedEndPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors, typename TPosition>
inline void
setClippedEndPosition(Gaps<TSequence, AnchorGaps<TGapAnchors> > & gaps, TPosition unclippedViewPosition)
{
    gaps.data_viewCutEnd = _unclippedLength(gaps) - unclippedViewPosition;
}

// ----------------------------------------------------------------------------
// Function clippedBeginPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors>
inline typename Position<Gaps<TSequence, AnchorGaps<TGapAnchors> > >::Type
clippedBeginPosition(Gaps<TSequence, AnchorGaps<TGapAnchors> > const & gaps)
{
    return gaps.data_viewCutBegin;
}

// ----------------------------------------------------------------------------
// Function clippedEndPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TGapAnchors>
inline typename Position<Gaps<TSequence, AnchorGaps<TGapAnchors> > >::Type
clippedEndPosition(Gaps<TSequence, AnchorGaps<TGapAnchors> > const & gaps)
{
    return _unclippedLength(gaps) - gaps.data_viewCutEnd;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ANCHOR_H_
