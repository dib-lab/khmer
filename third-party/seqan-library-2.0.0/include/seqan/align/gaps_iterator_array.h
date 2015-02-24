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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_ARRAY_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_ARRAY_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TGaps>
class Iter<TGaps, GapsIterator<ArrayGaps> >
{
public:
    // -----------------------------------------------------------------------
    // Internal Typedefs
    // -----------------------------------------------------------------------

    typedef typename TGaps::TArrayPos_        TArrayPos_;
    typedef typename TGaps::TArray_           TArray_;
    typedef typename Value<TArray_>::Type     TArrayValue_;
    typedef typename Position<TGaps>::Type    TGapsPos_;
    typedef typename Source<TGaps>::Type      TSource_;
    typedef typename Position<TSource_>::Type TSourcePos_;

    // -----------------------------------------------------------------------
    // Member Variables
    // -----------------------------------------------------------------------

    // The following index and position members must be mutable since the
    // insertion and deletion of gaps does not modify the iterator conceptually.

    // Pointer to the iterated container / gaps object.
    TGaps *     _container;
    // Index in the bucket array of the gaps object.
    mutable TArrayPos_  _bucketIndex;
    // Offset within the current bucket.
    mutable TArrayValue_  _bucketOffset;
    // Source position of the iterator.
    mutable TSourcePos_ _sourcePosition;
    // View position of the iterator.
    mutable TGapsPos_ _unclippedViewPosition;

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    // Default constructor.
    Iter() : _container(0), _bucketIndex(0), _bucketOffset(0), _sourcePosition(0), _unclippedViewPosition(0)
    {}

    // Copy constructor, required since we specify the one with complemented const below.
    Iter(Iter const & other) :
            _container(other._container), _bucketIndex(other._bucketIndex), _bucketOffset(other._bucketOffset),
            _sourcePosition(other._sourcePosition), _unclippedViewPosition(other._unclippedViewPosition)
    {}

    // Copy construtor for iterator -> const iterator conversion.
    // TODO(holtgrew): Think of something smarter, to restrict source types?
    template <typename TOtherIter>
    Iter(TOtherIter const & other) :
            _container(other._container), _bucketIndex(other._bucketIndex), _bucketOffset(other._bucketOffset),
            _sourcePosition(other._sourcePosition), _unclippedViewPosition(other._unclippedViewPosition)
    {}

    // Create at begin of gaps.
    Iter(TGaps & container_, Begin_ const &) :
            _container(&container_), _bucketIndex(0), _bucketOffset(0), _sourcePosition(0),
            _unclippedViewPosition(0)
    {
        if (_container->_array[0] == 0u)
            _bucketIndex = 1;
        // Go to beginning of clipping.
        goFurther(*this, container_._clippingBeginPos);
    }

    // Create at end of gaps.
    Iter(TGaps & container_, End_ const &) :
            _container(&container_), _bucketIndex(0), _bucketOffset(0), _sourcePosition(0),
            _unclippedViewPosition(0)
    {
        if (_container->_array[0] == 0u)
            _bucketIndex = 1;
        // Go to end of clipping position.
        goFurther(*this, container_._clippingEndPos);
    }

    // Create with position.
    // TODO(holtgrew): Chain constructor call to Begin_() here in C++11.
    template <typename TPos>
    Iter(TGaps & container_, TPos pos, Position_ const &) :
            _container(&container_), _bucketIndex(0), _bucketOffset(0), _sourcePosition(0),
            _unclippedViewPosition(0)
    {
        if (_container->_array[0] == 0u)
            _bucketIndex = 1;
        // pos is an unclipped view position, make it clipped.
        pos += container_._clippingBeginPos;
        // Initialized for begin position.  Now, advance.
        goFurther(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isGap()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
isGap(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    return !(it._bucketIndex % 2);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline typename Reference<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type
value(Iter<TGaps, GapsIterator<ArrayGaps> > & it)
{
    typedef typename Reference<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type TProxy;
    return TProxy(it);
}

template <typename TGaps>
inline typename Reference<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type
value(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    typedef typename Reference<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type TProxy;
    return TProxy(it);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Ideally, we would only have one here.

template <typename TGaps>
inline typename GetValue<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type
getValue(Iter<TGaps, GapsIterator<ArrayGaps> > & it)
{
    typedef typename Value<TGaps>::Type TAlphabet;
    if (isGap(it))
        return gapValue<TAlphabet>();
    else
        return value(source(container(it)), it._sourcePosition);
}

template <typename TGaps>
inline typename GetValue<Iter<TGaps, GapsIterator<ArrayGaps> > const>::Type
getValue(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    typedef typename Value<TGaps>::Type TAlphabet;
    if (isGap(it))
        return gapValue<TAlphabet>();
    else
        return value(source(container(it)), it._sourcePosition);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// Returns clipped view position of gaps iterator.

template <typename TGaps>
inline typename Position<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type
position(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    if (it._container == 0)
        return 0;

    typedef Iter<TGaps, GapsIterator<ArrayGaps> > TIter;
    typedef typename Position<TIter>::Type        TPosition;
    typedef typename TIter::TArrayPos_            TArrayPos;

    TPosition unclippedViewPosition = 0;
    for (TArrayPos i = 0; i < it._bucketIndex; ++i)
        unclippedViewPosition += it._container->_array[i];
    unclippedViewPosition += it._bucketOffset;

    // TODO(holtgrew): Simply return it._unclippedViewPosition?
    SEQAN_ASSERT_EQ(it._unclippedViewPosition, unclippedViewPosition);

    return unclippedViewPosition - clippedBeginPosition(*it._container);
}

// ----------------------------------------------------------------------------
// Function countGaps()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline typename Size<TGaps>::Type
countGaps(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    if (!isGap(it) || atEnd(it))
        return 0;  // Not on a gap or at end, no gap here.

    typedef typename Size<TGaps>::Type TSize;
    TSize result = it._container->_array[it._bucketIndex] - it._bucketOffset;
    // Check whether gaps reach behind the clipping and trim gaps for counting.
    if ((TSize)(it._unclippedViewPosition + result) > (TSize)it._container->_clippingEndPos)
        result = it._container->_clippingEndPos - it._unclippedViewPosition;
    return result;
}

// ----------------------------------------------------------------------------
// Function countCharacters()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline typename Size<TGaps>::Type
countCharacters(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    if (isGap(it) || atEnd(it))
        return 0;  // On a gap or at end, no characters here.

    typedef typename Size<TGaps>::Type TSize;
    TSize result = it._container->_array[it._bucketIndex] - it._bucketOffset;
    // Check whether gaps reach behind the clipping and trim gaps for counting.
    if ((TSize)(it._unclippedViewPosition + result) > (TSize)it._container->_clippingEndPos)
        result = it._container->_clippingEndPos - it._unclippedViewPosition;
    return result;
}

// ----------------------------------------------------------------------------
// Function goPrevious()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
goPrevious(Iter<TGaps, GapsIterator<ArrayGaps> > & it)
{
    typedef typename Position<TGaps>::Type TGapsPos;

    if (atBegin(it))  // Handle case of being at the beginning of the gaps.
        return false;

    if (it._bucketOffset > TGapsPos(0))
    {
        // Not at the beginning of a bucket.
        it._bucketOffset -= 1;
    }
    else
    {
        // At the beginning of a bucket.
        it._bucketIndex -= 1;
        SEQAN_ASSERT_GT(it._container->_array[it._bucketIndex], 0u);
        it._bucketOffset = it._container->_array[it._bucketIndex] - 1;
    }

    // Adjust source position.
    if (!isGap(it))
        it._sourcePosition -= 1;
    // Adjust clipped view position.
    it._unclippedViewPosition -= 1;

    return true;
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
goNext(Iter<TGaps, GapsIterator<ArrayGaps> > & it)
{
    if (atEnd(it))  // Handle case of being at the end of the gaps.
        return false;

    // Adjust source position.
    if (!isGap(it))
        it._sourcePosition += 1;
    // Adjust clipped view position.
    it._unclippedViewPosition += 1;

    if (it._bucketOffset + 1 != it._container->_array[it._bucketIndex])
    {
        // Not on last entry of a bucket.
        it._bucketOffset += 1;
    }
    else
    {
        // On last entry of a bucket.  If we are not in the last bucket then go
        // to next bucket.  Otherwise, go over the bucket, marks iterator-at-end.
        if (it._bucketIndex + 1 != length(it._container->_array))
        {
            // Go to next.
            it._bucketIndex += 1;
            if (it._bucketIndex > length(it._container->_array))
                SEQAN_ASSERT_GT(it._container->_array[it._bucketIndex], 0u);
            it._bucketOffset = 0;
        }
        else
        {
            // Go to end of bucket.
            it._bucketOffset += 1;
            SEQAN_ASSERT_EQ(it._bucketIndex + 1, length(it._container->_array));
            SEQAN_ASSERT_EQ(it._bucketOffset, back(it._container->_array));
        }
    }

    return true;
}


// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TDifference>
inline void
goFurther(Iter<TGaps, GapsIterator<ArrayGaps> > & it,
          TDifference delta)
{
    // TODO(holtgrew): Handle going backwards more efficiently.
    if (delta == TDifference(0))
        return;
    if (isNegative(delta))
    {
        typedef typename MakeSigned<TDifference>::Type TSignedDifference;
        for (; -static_cast<TSignedDifference>(delta); ++delta)
            goPrevious(it);
        return;
    }

    // Do nothing if we want are already at the end.
    if (atEnd(it))
        return;

    // Case: Going forward.

    // Get shortcut to new unclipped view position that we want to go to and
    // limit this to the clipping end pos of the gaps object.
    unsigned posEnd = it._unclippedViewPosition + delta;
    if (posEnd > static_cast<unsigned>(it._container->_clippingEndPos))
        posEnd = it._container->_clippingEndPos;

    // The variable counter is the number of view positions to go forward.
    for (unsigned counter = posEnd - it._unclippedViewPosition; counter > 0u;)
    {
        // Number of elements in bucket and number of remaining (behind offset in bucket).
        unsigned count = it._container->_array[it._bucketIndex];
        unsigned shift = count - it._bucketOffset;

        if (shift < counter)
        {
            it._unclippedViewPosition += shift;
            if (it._bucketIndex % 2)
                it._sourcePosition += shift;
            it._bucketOffset = 0;
            it._bucketIndex += 1;
            counter -= shift;
        }
        else if (shift == counter)
        {
            it._unclippedViewPosition += shift;
            if (it._bucketIndex % 2)
                it._sourcePosition += shift;

            // On last entry of a bucket.  If we are not in the last bucket then go to next bucket.
            // Otherwise, go over the bucket, marks iterator-at-end.
            if (it._bucketIndex + 1 != length(it._container->_array))
            {
                // Go to next.
                it._bucketIndex += 1;
                if (it._bucketIndex > length(it._container->_array))
                    SEQAN_ASSERT_GT(it._container->_array[it._bucketIndex], 0u);
                it._bucketOffset = 0;
            }
            else
            {
                // Go to end of bucket.
                it._bucketOffset += shift;
                SEQAN_ASSERT_EQ(it._bucketIndex + 1, length(it._container->_array));
                SEQAN_ASSERT_EQ(it._bucketOffset, back(it._container->_array));
            }
            counter = 0;
        }
        else  // shift > counter
        {
            it._unclippedViewPosition += counter;
            if (it._bucketIndex % 2)
                it._sourcePosition += counter;
            it._bucketOffset += counter;
            counter = 0;
        }
    }
}

// ----------------------------------------------------------------------------
// Function atBegin()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Non-const version is superflous :(
template <typename TGaps>
inline bool
atBegin(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    return it._unclippedViewPosition == it._container->_clippingBeginPos;
}

template <typename TGaps>
inline bool
atBegin(Iter<TGaps, GapsIterator<ArrayGaps> > & it)
{
    return it._unclippedViewPosition == it._container->_clippingBeginPos;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Non-const version is superflous :(
template <typename TGaps>
inline bool
atEnd(Iter<TGaps, GapsIterator<ArrayGaps> > const & it)
{
    return it._unclippedViewPosition == it._container->_clippingEndPos;
}

template <typename TGaps>
inline bool
atEnd(Iter<TGaps, GapsIterator<ArrayGaps> > & it)
{
    return it._unclippedViewPosition == it._container->_clippingEndPos;
}

// ----------------------------------------------------------------------------
// Function insertGaps()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TCount>
inline void
insertGaps(Iter<TGaps, GapsIterator<ArrayGaps> > const & it,
           TCount count)
{
    if (count == TCount(0))
        return;  // Do nothing!

    typedef typename TGaps::TArray_          TArray;
    typedef typename Position<TArray>::Type  TArrayPos;

    // Get shortcut to gaps.
    TGaps & gaps = *it._container;
    TArrayPos idx = it._bucketIndex;

    // Handle case of being at the start of a character bucket.
    if (idx % 2 && it._bucketOffset == 0)
    {
        idx -= 1;
        it._bucketIndex -= 1;
        it._bucketOffset = gaps._array[idx];
    }

    // Insert gaps, simple and fast if we are in a gaps bucket, a bit harder
    // otherwise.
    if (idx % 2)  // character bucket
    {
        if (gaps._array[idx] > it._bucketOffset)  // In the middle of the bucket.
        {
            TArray arr;
            resize(arr, 2, 0);
            arr[0] = count;
            arr[1] = gaps._array[idx] - it._bucketOffset;
            gaps._array[idx] = it._bucketOffset;
            insert(gaps._array, idx + 1, arr);

            // Update iterator.
            it._bucketIndex += 1;
            it._bucketOffset = 0;
        }
        else  // At the end of the bucket.
        {
            if (idx + 1 < length(gaps._array))  // Not at end of array.
            {
                gaps._array[idx + 1] += count;
            }
            else  // At end of array.
            {
                resize(gaps._array, length(gaps._array) + 2, 0);
                gaps._array[idx + 1] = count;
                gaps._array[idx + 2] = 0;
            }
        }
    }
    else  // gap bucket
    {
        gaps._array[idx] += count;
    }

    // Adjust clipping information.
    gaps._clippingEndPos += count;
}

// ----------------------------------------------------------------------------
// Function removeGaps()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TCount>
inline typename Size<TGaps>::Type
removeGaps(Iter<TGaps, GapsIterator<ArrayGaps> > const & it, TCount count)
{
    typedef typename TGaps::TArray_          TArray;
    typedef typename Position<TArray>::Type  TArrayPos;
    typedef typename Value<TArray>::Type     TArrayValue;
    typedef typename Source<TGaps>::Type     TSource;
    typedef typename Position<TSource>::Type TSeqPos;

    if (count == TCount(0))
        return 0;  // Do nothing!

    // Get shortcuts.
    TGaps & gaps = *it._container;
    TArrayPos idx = it._bucketIndex;
    TSeqPos offset = it._bucketOffset;

    // If we are inside a non-gap bucket then we cannot remove any gaps.
    if (idx % 2)
        return 0;

    // Otherwise, we can remove gaps right of the current position but not
    // more than there are.
    TSeqPos toRemove = count;
    if (toRemove > gaps._array[idx] - offset)
        toRemove = gaps._array[idx] - offset;
    gaps._array[idx] -= toRemove;
    // TODO(holtgrew): We have to decrement idx and adjust offset in case of merging!
    // In some cases, we remove the whole gap and merge the character buckets.
    if (gaps._array[idx] == TArrayValue(0))
    {
        // No merging for leading and trailing gap.
        if (idx == TArrayPos(0) || idx == TArrayPos(length(gaps._array) - 1))
        {
            gaps._array[idx - 1] += gaps._array[idx + 1];
            erase(gaps._array, idx, idx + 2);
        }
    }

    // Also update the right clipping position.
    gaps._clippingEndPos -= toRemove;

    // Finally, return number of removed gaps.
    return toRemove;
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
operator<(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs,
          Iter<TGaps, GapsIterator<ArrayGaps> > const & rhs)
{
    return lhs._bucketIndex < rhs._bucketIndex ||
            (lhs._bucketIndex == rhs._bucketIndex && lhs._bucketOffset < rhs._bucketOffset);
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
operator>(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs,
          Iter<TGaps, GapsIterator<ArrayGaps> > const & rhs)
{
    return lhs._bucketIndex > rhs._bucketIndex ||
            (lhs._bucketIndex == rhs._bucketIndex && lhs._bucketOffset > rhs._bucketOffset);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
operator<=(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs,
           Iter<TGaps, GapsIterator<ArrayGaps> > const & rhs)
{
    return !(lhs > rhs);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
operator>=(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs,
           Iter<TGaps, GapsIterator<ArrayGaps> > const & rhs)
{
    return !(lhs < rhs);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
operator==(Iter<TGaps, GapsIterator<ArrayGaps> > const & _lhs,
           Iter<TGaps, GapsIterator<ArrayGaps> > const & _rhs)
{
    return _lhs._container == _rhs._container &&
            _lhs._bucketIndex == _rhs._bucketIndex &&
            _lhs._bucketOffset == _rhs._bucketOffset;
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline bool
operator!=(Iter<TGaps, GapsIterator<ArrayGaps> > const & _lhs,
           Iter<TGaps, GapsIterator<ArrayGaps> > const & _rhs)
{
    return _lhs._container != _rhs._container ||
            _lhs._bucketIndex != _rhs._bucketIndex ||
            _lhs._bucketOffset != _rhs._bucketOffset;
}

// ----------------------------------------------------------------------------
// Function difference()
// ----------------------------------------------------------------------------

template <typename TGaps>
inline typename Difference<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type
difference(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs,
           Iter<TGaps, GapsIterator<ArrayGaps> > const & rhs)
{
    // TODO(holtgrew): Implementation could be more efficient.
    // We only need to solve the case lhs < rhs.
    if (lhs > rhs)
        return -difference(rhs, lhs);
    if (lhs == rhs)
        return 0;
    SEQAN_ASSERT(lhs < rhs);  // Makes code below simpler.

    typedef Iter<TGaps, GapsIterator<ArrayGaps> > TIter;
    typedef typename Difference<TIter>::Type      TDifference;
    TDifference d = 0;
    for (TIter it = lhs; it != rhs; ++it)
        ++d;

    return -d;
}

// ----------------------------------------------------------------------------
// Function operator-()                                            [difference]
// ----------------------------------------------------------------------------

template <typename TGaps>
inline typename Difference<Iter<TGaps, GapsIterator<ArrayGaps> > >::Type
operator-(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs,
          Iter<TGaps, GapsIterator<ArrayGaps> > const & rhs)
{
    return difference(lhs, rhs);
}

// ----------------------------------------------------------------------------
// Function operator-()                                         [copy movement]
// ----------------------------------------------------------------------------

template <typename TGaps, typename TDifference>
inline Iter<TGaps, GapsIterator<ArrayGaps> >
operator-(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs, TDifference d)
{
    Iter<TGaps, GapsIterator<ArrayGaps> > result = lhs;
    result -= d;
    return result;
}

// ----------------------------------------------------------------------------
// Function operator+()                                         [copy movement]
// ----------------------------------------------------------------------------

template <typename TGaps, typename TDifference>
inline Iter<TGaps, GapsIterator<ArrayGaps> >
operator+(Iter<TGaps, GapsIterator<ArrayGaps> > const & lhs, TDifference d)
{
    Iter<TGaps, GapsIterator<ArrayGaps> > result = lhs;
    result += d;
    return result;
}

}  // namespace seqan

#endif  // SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_ARRAY_H_
