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

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_ANCHOR_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_ANCHOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Iterator for AnchorGaps objects.

template <typename TGaps_, typename TGapAnchors_> //Gaps<TSource_, AnchorGaps<TGapAnchors_> >
class Iter<TGaps_, GapsIterator<AnchorGaps<TGapAnchors_> > >
{
public:
    typedef TGaps_                                            TGaps;
    typedef typename Source<TGaps>::Type                    TSource;
    typedef TGapAnchors_                                    TGapAnchors;

    // TODO(holtgrew): Why is the following commented out?
//    typedef typename Value<TGapAnchors>::Type                TGapAnchor;
    typedef typename Size<typename Value<TGapAnchors>::Type>::Type          TGapAnchorSize_;
    typedef GapAnchor<typename MakeSigned_<TGapAnchorSize_>::Type>          TGapAnchor;
    typedef typename MakeSigned<typename Position<TGapAnchor>::Type>::Type  TGapPos;
    typedef typename Iterator<TGapAnchors, Standard>::Type                  TAnchorIter;

    TGaps *                    data_container;                            //the gaps object
    TGapPos                 seqLength;
    mutable TGapAnchor        current;
    mutable TGapAnchor        prevAnchor;
    mutable TGapAnchor        nextAnchor;
    mutable TGapAnchor        viewBegin;
    mutable TGapAnchor        viewEnd;
    mutable int                anchorIdx;

public:
    Iter()
    {
SEQAN_CHECKPOINT
        data_container = NULL;
        seqLength = 0;
    }
/*    Iter(Iter const & other_):
        data_container(other_.data_container),
        seqLength(other_.seqLength),
        current(other_.current),
        prevAnchor(other_.prevAnchor),
        nextAnchor(other_.nextAnchor),
        anchorIdx(other_.anchorIdx)
    {
SEQAN_CHECKPOINT
    }
*/    Iter(TGaps & container_):
        data_container(&container_)
    {
SEQAN_CHECKPOINT
        _assignSourceLength(seqLength, container_);
        _goToGapAnchorIterator(*this, data_container->data_viewCutBegin + data_container->data_cutBegin);
        viewBegin = current;
        viewEnd.gapPos   = _unclippedLength(*data_container) + data_container->data_cutBegin - data_container->data_viewCutEnd;
        viewEnd.seqPos = positionGapToSeq(*data_container, viewEnd.gapPos);
    }
    Iter(TGaps & container_, TGapPos clippedViewPosition):
        data_container(&container_)
    {
SEQAN_CHECKPOINT
        _assignSourceLength(seqLength, container_);
        _goToGapAnchorIterator(*this, clippedViewPosition + data_container->data_viewCutBegin + data_container->data_cutBegin);
        viewBegin.gapPos = data_container->data_viewCutBegin + data_container->data_cutBegin;
        viewEnd.gapPos   = _unclippedLength(*data_container) + data_container->data_cutBegin - data_container->data_viewCutEnd;
        viewBegin.seqPos = positionGapToSeq(*data_container, viewBegin.gapPos);
        viewEnd.seqPos   = positionGapToSeq(*data_container, viewEnd.gapPos);
    }
    ~Iter()
    {
SEQAN_CHECKPOINT
    }

    Iter const & operator = (Iter const & other_)
    {
SEQAN_CHECKPOINT
        data_container = other_.data_container;
        seqLength = other_.seqLength;
        current = other_.current;
        prevAnchor = other_.prevAnchor;
        nextAnchor = other_.nextAnchor;
        anchorIdx = other_.anchorIdx;
        viewBegin = other_.viewBegin;
        viewEnd = other_.viewEnd;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Can go if data_container were _container, dupe in _base.h

template <typename TGaps, typename TGapsArray>
inline TGaps &
container(Iter<TGaps, GapsIterator<AnchorGaps<TGapsArray> > > & me)
{
    return *me.data_container;
}

template <typename TGaps, typename TGapsArray>
inline TGaps &
container(Iter<TGaps, GapsIterator<AnchorGaps<TGapsArray> > > const & me)
{
    return *me.data_container;
}

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Source<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type
source(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
    return begin(source(*me.data_container), Rooted()) + me.current.seqPos;
}

template <typename TGaps, typename TGapAnchors>
inline typename Source<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
source(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    return begin(source(*me.data_container), Rooted()) + me.current.seqPos;
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename GetValue< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
getValue(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
    typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type TValue;
    if (isGap(me)) return gapValue<TValue>();
    else if (isUnknown(me)) return unknownValue<TValue>();
    else return getValue(source(me));
}

template <typename TGaps, typename TGapAnchors>
inline typename GetValue< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type
getValue(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    typedef typename Value<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type TValue;
    if (isGap(me)) return gapValue<TValue>();
    else if (isUnknown(me)) return unknownValue<TValue>();
    else return getValue(source(me));
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Reference< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
value(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & it)
{
    typedef typename Reference<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type TProxy;
    return TProxy(it);
}

template <typename TGaps, typename TGapAnchors>
inline typename Reference< Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type
value(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & it)
{
    typedef typename Reference<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const>::Type TProxy;
    return TProxy(it);
}

// ----------------------------------------------------------------------------
// Function isGap()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
isGap(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    return me.current.seqPos == me.nextAnchor.seqPos;
}

// ----------------------------------------------------------------------------
// Function isUnknown()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
isUnknown(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    int len;
    _assignSourceLength(len, *me.data_container);
    return me.current.seqPos < 0 || me.current.seqPos >= len;
}

// ----------------------------------------------------------------------------
// Function isClipped()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
isClipped(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    return me.current.gapPos == me.nextAnchor.gapPos;
}

// ----------------------------------------------------------------------------
// Function countGaps()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Size<TGaps>::Type
countGaps(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    if (!isGap(me))
        return 0;
    if (me.nextAnchor.gapPos > me.viewEnd.gapPos)
        return me.viewEnd.gapPos - me.current.gapPos;
    return me.nextAnchor.gapPos - me.current.gapPos;
}

// ----------------------------------------------------------------------------
// Function countCharacters()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Size<TGaps>::Type
countCharacters(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
    if (isGap(me))
        return 0;
    if (me.nextAnchor.seqPos > me.viewEnd.seqPos)
        return me.viewEnd.seqPos - me.current.seqPos;
    return me.nextAnchor.seqPos - me.current.seqPos;
}

// ----------------------------------------------------------------------------
// Function blockLength()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Size<TGaps>::Type
blockLength(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
    if (isGap(me))
        return countGaps(me);
    else
        return countCharacters(me);
}

// ----------------------------------------------------------------------------
// Function atBegin()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
atBegin(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//    return me.current.seqPos == 0 && me.current.gapPos == 0;
    return me.current <= me.viewBegin;
}

template <typename TGaps, typename TGapAnchors>
inline bool
atBegin(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//    return me.current.seqPos == 0 && me.current.gapPos == 0;
    return me.current <= me.viewBegin;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
atEnd(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
//    return me.current == me.nextAnchor;
    return me.current >= me.viewEnd;
}

template <typename TGaps, typename TGapAnchors>
inline bool
atEnd(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me)
{
//    return me.current == me.nextAnchor;
    return me.current >= me.viewEnd;
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
operator == (
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
    return left.current == right.current;
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
operator != (
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
    return left.current != right.current;
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
operator < (
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
    return left.current < right.current;
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
operator<=(
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
    return !(left.current > right.current);
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
operator > (
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & left,
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & right)
{
    return left.current > right.current;
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline bool
operator>=(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & lhs,
           Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & rhs)
{
    return !(lhs < rhs);
}

// ----------------------------------------------------------------------------
// Function insertGaps()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors, typename TCount>
inline void
insertGaps(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & me,
           TCount size)
{
    TGapAnchors & anchors = _dataAnchors(*me.data_container);
    typedef typename Iterator<TGapAnchors, Standard>::Type TIter;

    if (size <= 0) return;

    // insert a new anchor
    if (!isGap(me))
    {
        if (me.prevAnchor.gapPos == me.current.gapPos)
        {
            me.nextAnchor = me.prevAnchor;
            _getAnchor(me.prevAnchor, *me.data_container, --me.anchorIdx);
        }
        else
        {
            me.nextAnchor = me.current;
            insertValue(anchors, me.anchorIdx, me.nextAnchor, Generous());
        }
    }
    else
    {
        if (me.anchorIdx >= (int)length(anchors))
        {
            // add gap after the sequence and in (or at the right boundary of) the view
            if (me.current.gapPos <= me.viewEnd.gapPos)
            {
                container(me).data_cutEnd -= size;
                me.viewEnd.gapPos += size;
            }
            return;
        }
        if (empty(anchors))
            appendValue(anchors, me.nextAnchor, Generous());
    }
    if (me.anchorIdx < (int)length(anchors))
    {
        if (me.anchorIdx >= 0)
        {
            me.nextAnchor.gapPos += size;
            TIter it = begin(anchors, Standard());
            TIter itEnd = end(anchors, Standard());
            if (me.anchorIdx >= 0)
                it += me.anchorIdx;
            for (; it != itEnd; ++it)
                (*it).gapPos += size;
        }
        else
            // add gap before the sequence and in (or at the left boundary of) the view
            if (me.current.gapPos >= me.viewBegin.gapPos)
            {
                container(me).data_cutBegin -= size;
                me.viewBegin.gapPos -= size;
                me.current.gapPos -= size;
                return;
            }
    }
    if (me.current.gapPos <= me.viewEnd.gapPos)
        me.viewEnd.gapPos += size;

/*
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > it2 = begin(*me.data_container) + me.current.gapPos;
    if (me.current != it2.current || me.prevAnchor != it2.prevAnchor || me.nextAnchor != it2.nextAnchor || me.anchorIdx != it2.anchorIdx)
        std::cout<<"*";
*/
}

// ----------------------------------------------------------------------------
// Function removeGaps()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors, typename TCount>
inline typename Size<TGaps>::Type
removeGaps(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & it,
           TCount size_)
{
    TGapAnchors & anchors = _dataAnchors(*it.data_container);
    typedef typename Iterator<TGapAnchors, Standard>::Type TAnchorsIter;

    typedef Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > TIter;
    typedef typename TIter::TGapAnchor TGapAnchor;
    // typedef typename Value<TGapAnchors>::Type   TGapAnchor;
    typedef typename Position<TGapAnchor>::Type TPos;

    if (size_ <= 0 || !isGap(it))
        return 0;
    TPos size = size_;

    // static_cast<TGapAnchor>(Nothing());
    // static_cast<TPos>(Nothing());
    if (it.current.gapPos + size > it.nextAnchor.gapPos)
        size = it.nextAnchor.gapPos - it.current.gapPos;

    if (it.prevAnchor.gapPos + it.current.seqPos == it.current.gapPos + it.prevAnchor.seqPos &&
        it.current.gapPos + size == it.nextAnchor.gapPos)
    {
        // remove the gap
        if (it.anchorIdx < (int)length(anchors))
            erase(anchors, it.anchorIdx);
        _getAnchor(it.nextAnchor, *it.data_container, it.anchorIdx);
    }

    // shift anchors
    if (it.anchorIdx < (int)length(anchors))
    {
        if (it.anchorIdx >= 0)
        {
            it.nextAnchor.gapPos -= size;
            TAnchorsIter itA = begin(anchors, Standard());
            TAnchorsIter itAEnd = end(anchors, Standard());
            if (it.anchorIdx >= 0)
                itA += it.anchorIdx;
            for (; itA != itAEnd; ++itA)
                (*itA).gapPos -= size;
        }
        else
            // remove gap before the sequence and in (or at the left boundary of) the view
            if (it.current.gapPos >= it.viewBegin.gapPos)
            {
                // assure that we don't remove more gaps than available
                if (size > it.nextAnchor.gapPos - it.current.gapPos)
                    size = it.nextAnchor.gapPos - it.current.gapPos;
                container(it).data_cutBegin += size;
                it.viewBegin.gapPos += size;
                it.current.gapPos += size;
                return size;
            }
    }
    else
    {
        if (it.current.gapPos <= it.viewEnd.gapPos)
            container(it).data_cutEnd += size;
    }
    if (it.current.gapPos <= it.viewEnd.gapPos)
        it.viewEnd.gapPos -= size;

    return size;
/*
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > it2 = begin(*me.data_container) + me.current.gapPos;
    if (me.current != it2.current || me.prevAnchor != it2.prevAnchor || me.nextAnchor != it2.nextAnchor || me.anchorIdx != it2.anchorIdx)
        std::cout<<"*";
*/
}

// ----------------------------------------------------------------------------
// Helper Function _goNextGapAnchorIterator()
// ----------------------------------------------------------------------------

template <typename T>
inline void
_goNextGapAnchorIterator(T & me)
{
    if (me.current.gapPos < me.nextAnchor.gapPos)
    {
        ++me.current.gapPos;
        if (me.current.seqPos < me.nextAnchor.seqPos)
            ++me.current.seqPos;
    }
    while (me.current.gapPos == me.nextAnchor.gapPos)
    {
        me.current = me.prevAnchor = me.nextAnchor;
        _getAnchor(me.nextAnchor, *me.data_container, ++me.anchorIdx + 1);
    }
}

// ----------------------------------------------------------------------------
// Helper Function _goPreviousGapAnchorIterator()
// ----------------------------------------------------------------------------

template <typename T>
inline void
_goPreviousGapAnchorIterator(T & me)
{
    while (me.current.gapPos == me.prevAnchor.gapPos)
    {
        me.current = me.nextAnchor = me.prevAnchor;
        _getAnchor(me.prevAnchor, *me.data_container, --me.anchorIdx);
    }
    --me.current.gapPos;
    if (me.nextAnchor.seqPos - me.prevAnchor.seqPos > me.current.gapPos - me.prevAnchor.gapPos)
        me.current.seqPos = me.prevAnchor.seqPos + (me.current.gapPos - me.prevAnchor.gapPos);
    else
        me.current.seqPos = me.nextAnchor.seqPos;
}

// ----------------------------------------------------------------------------
// Helper Function _goToGapAnchorIterator()
// ----------------------------------------------------------------------------

template <typename T, typename TPos>
inline void
_goToGapAnchorIterator(T & me, TPos pos)
{
    typedef typename T::TGapAnchors                 TGapAnchors;
    typedef typename Value<TGapAnchors>::Type       TGapAnchor;
    typedef typename Position<TGapAnchor>::Type     TAnchorPos;
    typedef typename MakeSigned<TAnchorPos>::Type   TAnchorSPos;

    if (isNegative(pos))
        me.anchorIdx = -1;
    else
    {
        TGapAnchors const & anchors = _dataAnchors(*me.data_container);
        if (!empty(anchors))
        {
            me.anchorIdx = upperBoundGapAnchor(anchors, pos, SortGapPos()) - begin(anchors, Standard());
            if (me.anchorIdx < (int)length(anchors))
                if (anchors[me.anchorIdx].gapPos == (TAnchorPos)pos && anchors[me.anchorIdx].seqPos != (TAnchorPos)me.seqLength)
                    ++me.anchorIdx;
        }
        else
        {
            me.anchorIdx = ((TAnchorSPos)pos < me.seqLength)? 0: 1;
        }
    }
    _getAnchor(me.prevAnchor, *me.data_container, me.anchorIdx);
    _getAnchor(me.nextAnchor, *me.data_container, me.anchorIdx + 1);

    me.current.gapPos = pos;
    if (me.nextAnchor.seqPos - me.prevAnchor.seqPos > (int)pos - me.prevAnchor.gapPos)
        me.current.seqPos = me.prevAnchor.seqPos + ((int)pos - me.prevAnchor.gapPos);
    else
        me.current.seqPos = me.nextAnchor.seqPos;
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline void
goNext(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
    _goNextGapAnchorIterator(me);
}

// ----------------------------------------------------------------------------
// Function goPrevious()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline void
goPrevious(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me)
{
    _goPreviousGapAnchorIterator(me);
}

// ----------------------------------------------------------------------------
// Function goFurther()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors, typename TSize>
inline void
goFurther(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > & me, TSize steps)
{
    _goToGapAnchorIterator(me, me.current.gapPos + steps);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// Returns clipped view position of gaps iterator.

template <typename TGaps, typename TGapAnchors>
inline typename Position<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
position(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & it)
{
    return it.current.gapPos - it.viewBegin.gapPos;
}

// ----------------------------------------------------------------------------
// Function difference()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Difference<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
difference(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & lhs,
           Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & rhs)
{
    return lhs.current.gapPos - rhs.current.gapPos;
}

// ----------------------------------------------------------------------------
// Function operator-()                                            [difference]
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors>
inline typename Difference<Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > >::Type
operator-(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & lhs,
          Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & rhs)
{
    return difference(lhs, rhs);
}

// ----------------------------------------------------------------------------
// Function operator-()                                         [copy movement]
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors, typename TDifference>
inline Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > >
operator-(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & lhs, TDifference d)
{
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > result = lhs;
    result -= d;
    return result;
}

// ----------------------------------------------------------------------------
// Function operator+()                                         [copy movement]
// ----------------------------------------------------------------------------

template <typename TGaps, typename TGapAnchors, typename TDifference>
inline Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > >
operator+(Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > const & lhs, TDifference d)
{
    Iter<TGaps, GapsIterator<AnchorGaps<TGapAnchors> > > result = lhs;
    result += d;
    return result;
}


}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_ANCHOR_H_
