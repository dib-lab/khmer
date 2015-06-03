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

// SEQAN_NO_GENERATED_FORWARDS

// TODO(holtgrew): Currently, operations are a function of the whole gap count, could be of clipped region only.
// TODO(holtgrew): Problem with the gap value, getValue(), value().

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ARRAY_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ARRAY_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// Internally used tag for creating iterators at the begin of containers.
struct Begin__;
typedef Tag<Begin__> Begin_;

// Internally used tag for creating iterators at the end of containers.
struct End__;
typedef Tag<End__> End_;

// Internally used tag for creating iterators inside of containers.
struct Position__;
typedef Tag<Position__> Position_;

struct ArrayGaps_;
typedef Tag<ArrayGaps_> ArrayGaps;

template <typename TSequence> class Gaps<TSequence, ArrayGaps>;

template <typename TSequence>
inline void _reinitArrayGaps(Gaps<TSequence, ArrayGaps> & gaps);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct ArrayGaps_;
typedef Tag<ArrayGaps_> ArrayGaps;

/*!
 * @class ArrayGaps
 * @headerfile <seqan/align.h>
 * @extends Gaps
 * @brief Stores length of gap- and non-gap runs in an array.
 *
 * @signature template <typename TSequence>
 *            class Gaps<TSequence, ArrayGaps>
 *
 * @tparam TSequence The type of the underling sequence.
 */

/*!
 * @fn ArrayGaps::Gaps
 * @headerfile <seqan/align.h>
 * @brief Constructor.
 *
 * @signature Gaps::Gaps([other]);
 * @signature Gaps::Gaps(seq);
 *
 * @param[in] other Other Gaps object to copy from.
 * @param[in] seq   Sequence concept to construct the gaps for.
 */

template <typename TSequence>
class Gaps<TSequence, ArrayGaps>
{
public:
    // -----------------------------------------------------------------------
    // Internal Typedefs
    // -----------------------------------------------------------------------

    typedef typename Size<Gaps>::Type          TSize_;
    typedef typename Size<TSequence>::Type     TSequenceSize_;
    typedef typename Position<Gaps>::Type      TPosition_;
    typedef typename Position<TSequence>::Type TSequencePosition_;
    typedef typename Value<Gaps>::Type         TValue_;

    typedef String<TSequenceSize_>             TArray_;
    typedef typename Position<TArray_>::Type   TArrayPos_;

    // -----------------------------------------------------------------------
    // Member Variables
    // -----------------------------------------------------------------------

    // Holder of the underlying sequence.
    Holder<TSequence> _source;

    // The array with the alternating gap/source char counts.
    TArray_ _array;

    // Begin and end position in the source.
    TSequencePosition_ _sourceBeginPos, _sourceEndPos;
    // Begin and end position in the view.
    TPosition_ _clippingBeginPos, _clippingEndPos;
    // TODO(holtgrew): The following is a possible optimization.
    // // Index of clipping begin and end in the _array seq/gap char count array.
    // // This identifies a slice of the view.
    // TArrayPos_ _clippingBeginIdx, _clippingEndIdx;
    // // Offset within the slice.
    // TSequenceSize_ _clippingBeginOffset, _clippingEndOffset;

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    Gaps() : _sourceBeginPos(0), _sourceEndPos(0), _clippingBeginPos(0), _clippingEndPos(0)//,
             // _clippingBeginIdx(0), _clippingEndIdx(0), _clippingBeginOffset(0), _clippingEndOffset(0)
    {}

    explicit
    Gaps(TSequence & seq) :
            _source(seq), _sourceBeginPos(0), _sourceEndPos(length(seq)),
            _clippingBeginPos(0), _clippingEndPos(length(seq))//,
            // _clippingBeginIdx(0), _clippingEndIdx(0), _clippingBeginOffset(0),
            // _clippingEndOffset(0)
    {
        // Initialize array gaps object for ungapped sequence.
        _reinitArrayGaps(*this);
    }

    Gaps(Gaps const & other) :
        _source(other._source), _array(other._array), _sourceBeginPos(other._sourceBeginPos),
        _sourceEndPos(other._sourceEndPos), _clippingBeginPos(other._clippingBeginPos),
        _clippingEndPos(other._clippingEndPos)
    {}

    // -----------------------------------------------------------------------
    // Array Subscript Operator
    // -----------------------------------------------------------------------

    inline Gaps &
    operator=(Gaps const & other)
    {
        setValue(_source, source(other));
        _array = other._array;
        _sourceBeginPos = other._sourceBeginPos;
        _sourceEndPos = other._sourceEndPos;
        _clippingBeginPos = other._clippingBeginPos;
        _clippingEndPos = other._clippingEndPos;
        return *this;
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

// ----------------------------------------------------------------------------
// Function swap()
// ----------------------------------------------------------------------------

template <typename TSequence>
void swap(Gaps<TSequence, ArrayGaps> & lhs, Gaps<TSequence, ArrayGaps> & rhs)
{
    swap(lhs._source, rhs._source);
    swap(lhs._array, rhs._array);

    std::swap(lhs._sourceBeginPos, rhs._sourceBeginPos);
    std::swap(lhs._sourceEndPos, rhs._sourceEndPos);
    std::swap(lhs._clippingBeginPos, rhs._clippingBeginPos);
    std::swap(lhs._clippingEndPos, rhs._clippingEndPos);
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function detach()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove? Only used by module blast.

template <typename TSequence>
void detach(Gaps<TSequence, ArrayGaps> & gaps)
{
    detach(gaps._source);
}

// ----------------------------------------------------------------------------
// Function _setLength()
// ----------------------------------------------------------------------------

// Set the length, only use if TSequence is Nothing.

template <typename TSequence, typename TSize>
inline void _setLength(Gaps<TSequence, ArrayGaps> & gaps, TSize newLen)
{
    // Reset array.
    resize(gaps._array, 3);
    gaps._array[0] = 0;
    gaps._array[1] = newLen;
    gaps._array[2] = 0;
    // Reset clipping information.
    gaps._clippingBeginPos = 0;
    gaps._clippingEndPos = newLen;
    gaps._sourceBeginPos = 0;
    gaps._sourceEndPos = gaps._clippingEndPos;
    // gaps._clippingBeginIdx = 1;
    // gaps._clippingBeginOffset = 0;
    // gaps._clippingEndIdx = 1;
    // gaps._clippingEndOffset = value(gaps._source)[1];
}

// ----------------------------------------------------------------------------
// Helper Function _reinitArrayGaps()
// ----------------------------------------------------------------------------

// Reset the array gaps DS such that represents the ungapped sequence.

template <typename TSequence>
inline void _reinitArrayGaps(Gaps<TSequence, ArrayGaps> & gaps)
{
    _setLength(gaps, length(value(gaps._source)));
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We'd rather have "TTag const &" here.
template <typename TSequence, typename TTag>
inline typename Iterator<Gaps<TSequence, ArrayGaps> >::Type
begin(Gaps<TSequence, ArrayGaps> & gaps, Tag<TTag> const /*tag*/)
{
    typedef typename Iterator<Gaps<TSequence, ArrayGaps> >::Type TIter;
    return TIter(gaps, Begin_());
}

// TODO(holtgrew): We'd rather have "TTag const &" here.
template <typename TSequence, typename TTag>
inline typename Iterator<Gaps<TSequence, ArrayGaps> const>::Type
begin(Gaps<TSequence, ArrayGaps> const & gaps, Tag<TTag> const /*tag*/)
{
    typedef typename Iterator<Gaps<TSequence, ArrayGaps> const>::Type TIter;
    return TIter(gaps, Begin_());
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We'd rather have "TTag const &" here.
template <typename TSequence, typename TTag>
inline typename Iterator<Gaps<TSequence, ArrayGaps> >::Type
end(Gaps<TSequence, ArrayGaps> & gaps, Tag<TTag> const /*tag*/)
{
    typedef typename Iterator<Gaps<TSequence, ArrayGaps> >::Type TIter;
    return TIter(gaps, End_());
}

// TODO(holtgrew): We'd rather have "TTag const &" here.
template <typename TSequence, typename TTag>
inline typename Iterator<Gaps<TSequence, ArrayGaps> const>::Type
end(Gaps<TSequence, ArrayGaps> const & gaps, Tag<TTag> const /*tag*/)
{
    typedef typename Iterator<Gaps<TSequence, ArrayGaps> const>::Type TIter;
    return TIter(gaps, End_());
}

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We'd rather have "TTag const &" here.
template <typename TSequence, typename TTag, typename TPosition>
inline typename Iterator<Gaps<TSequence, ArrayGaps> >::Type
iter(Gaps<TSequence, ArrayGaps> & gaps, TPosition pos, Tag<TTag> const /*tag*/)
{
    typedef typename Iterator<Gaps<TSequence, ArrayGaps> >::Type TIter;
    return TIter(gaps, pos, Position_());
}

// TODO(holtgrew): We'd rather have "TTag const &" here.
template <typename TSequence, typename TTag, typename TPosition>
inline typename Iterator<Gaps<TSequence, ArrayGaps> const>::Type
iter(Gaps<TSequence, ArrayGaps> const & gaps, TPosition pos, Tag<TTag> const /*tag*/)
{
    typedef typename Iterator<Gaps<TSequence, ArrayGaps> const>::Type TIter;
    return TIter(gaps, pos, Position_());
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline typename Size<Gaps<TSequence, ArrayGaps> >::Type
length(Gaps<TSequence, ArrayGaps> const & gaps)
{
    SEQAN_ASSERT_GEQ(gaps._clippingEndPos, gaps._clippingBeginPos);
    return gaps._clippingEndPos - gaps._clippingBeginPos;
}

// ----------------------------------------------------------------------------
// Function unclippedLength()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline typename Size<Gaps<TSequence, ArrayGaps> >::Type
unclippedLength(Gaps<TSequence, ArrayGaps> const & gaps)
{
    typedef typename Size<Gaps<TSequence, ArrayGaps> >::Type TSize;

    TSize result = 0;
    for (unsigned i = 0; i < length(gaps._array); ++i)
        result += gaps._array[i];

    return result;
}

// ----------------------------------------------------------------------------
// Function createSource()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove? Switch to Hosted Type Interface?

template <typename TSequence>
inline void createSource(Gaps<TSequence, ArrayGaps> & gaps)
{
    create(gaps._source);
}

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline typename Source<Gaps<TSequence, ArrayGaps> const>::Type &
source(Gaps<TSequence, ArrayGaps> const & gaps)
{
    return value(gaps._source);
}

template <typename TSequence>
inline typename Source<Gaps<TSequence, ArrayGaps> >::Type &
source(Gaps<TSequence, ArrayGaps> & gaps)
{
    return value(gaps._source);
}

// ----------------------------------------------------------------------------
// Function setSource()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Test with clippings, also for AnchorGaps.

template <typename TSequence>
inline void
setSource(Gaps<TSequence, ArrayGaps> & gaps, TSequence & source)
{
    setValue(gaps._source, source);
    _reinitArrayGaps(gaps);
}

template <typename TSequence>
inline void
setSource(Gaps<TSequence const, ArrayGaps> & gaps, TSequence & source)
{
    setValue(gaps._source, source);
    _reinitArrayGaps(gaps);
}

// ----------------------------------------------------------------------------
// Function assignSource()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSequence2>
inline void
assignSource(Gaps<TSequence, ArrayGaps> & gaps, TSequence2 const & source)
{
    value(gaps._source) = source;
    _reinitArrayGaps(gaps);
}

// ----------------------------------------------------------------------------
// Function toSourcePosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline typename Position<TSequence>::Type
toSourcePosition(Gaps<TSequence, ArrayGaps> const & gaps, TPosition clippedViewPos)
{
    typedef Gaps<TSequence, ArrayGaps>         TGaps;
    typedef typename Position<TGaps>::Type     TGapsPos;
    typedef typename TGaps::TArrayPos_         TArrayPos;
    typedef typename Position<TSequence>::Type TSourcePos;

    // Translate from clipped view position to unclipped view position.
    TGapsPos unclippedViewPos = clippedViewPos + clippedBeginPosition(gaps);

    // Get index i of the according bucket and offset within bucket.
    TSourcePos result = 0;
    TArrayPos i = 0;
    TSourcePos const iEnd = length(gaps._array);
    for (TSourcePos counter = unclippedViewPos; counter > TGapsPos(0) && i < iEnd;)
    {
        if (counter > gaps._array[i])
        {
            if (i % 2)  // character bucket
                result += gaps._array[i];
            counter -= gaps._array[i];
            i += 1;
        }
        else if (counter <= gaps._array[i])
        {
            if (i % 2)  // character bucket
            {
                result += counter;
            }
            counter = 0;
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
// Function toViewPosition()
// ----------------------------------------------------------------------------

// Parameter rightOfGaps moves to the right end of gaps if the character at sourcePosition is followed by a gap in the
// view.
template <typename TSequence, typename TPosition>
inline typename Position<Gaps<TSequence, ArrayGaps> >::Type
toViewPosition(Gaps<TSequence, ArrayGaps> const & gaps, TPosition sourcePosition, bool rightOfGaps = true)
{
    typedef Gaps<TSequence, ArrayGaps>     TGaps;
    typedef typename Position<TGaps>::Type TGapsPosition;
    typedef typename TGaps::TArray_        TArray;
    typedef typename TGaps::TArrayPos_     TArrayPos;
    typedef typename Value<TArray>::Type   TArrayValue;

    if (sourcePosition == TPosition(0))
        return gaps._array[0] - clippedBeginPosition(gaps);

    // First, convert to unclipped source position.
    TGapsPosition unclippedViewPosition = 0;
    TArrayPos i = 0;
    for (TArrayValue counter = sourcePosition; counter > TArrayValue(0); ++i)
    {
        if (i % 2 /*== 1*/)  // sequence bucket
        {
            if (counter > gaps._array[i])
            {
                unclippedViewPosition += gaps._array[i];
                counter -= gaps._array[i];
            }
            else if (counter < gaps._array[i])
            {
                unclippedViewPosition += counter;
                counter = 0;
            }
            else  // counter == gaps._array[i]
            {
                unclippedViewPosition += counter;
                if (rightOfGaps && i + 2 < length(gaps._array))
                    unclippedViewPosition += gaps._array[i + 1];
                counter = 0;
            }
        }
        else  // gaps bucket
        {
            unclippedViewPosition += gaps._array[i];
        }
    }

    // Return after clipping.
    return unclippedViewPosition - clippedBeginPosition(gaps);
}

// ----------------------------------------------------------------------------
// Function insertGaps()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition, typename TCount>
inline void
insertGaps(Gaps<TSequence, ArrayGaps> & gaps, TPosition clippedViewPos, TCount count)
{
    typedef Gaps<TSequence, ArrayGaps>     TGaps;
    typedef typename Position<TGaps>::Type TGapsPosition;
    typedef typename TGaps::TArray_        TArray;
    typedef typename TGaps::TArrayPos_     TArrayPos;
    typedef typename Position<TSequence>::Type TSeqPos;

    // Translate from clipped view position to unclipped view position.
    TGapsPosition unclippedViewPos = clippedViewPos + clippedBeginPosition(gaps);

    // Get index i of the according bucket and offset within bucket.
    TArrayPos i = 0;
    TSeqPos offset = 0;
    for (TSeqPos counter = unclippedViewPos; counter > 0;)
    {
        SEQAN_ASSERT_LT(i, length(gaps._array));
        if (counter > gaps._array[i])
        {
            counter -= gaps._array[i];
            i += 1;
        }
        else
        {
            offset = counter;
            counter = 0;
        }
    }

    SEQAN_ASSERT_GEQ(gaps._array[i], offset);

    // Insert gaps, simple and fast if we are in a gaps bucket, a bit harder
    // otherwise.
    if (i % 2)  // character bucket
    {
        if (gaps._array[i] > offset)  // In the middle of the bucket.
        {
            TArray arr;
            resize(arr, 2, 0);
            arr[0] = count;
            arr[1] = gaps._array[i] - offset;
            gaps._array[i] = offset;
            insert(gaps._array, i + 1, arr);
        }
        else  // At the end of the bucket.
        {
            if (i + 1 < length(gaps._array))  // Not at end of array.
            {
                gaps._array[i + 1] += count;
            }
            else  // At end of array.
            {
                resize(gaps._array, length(gaps._array) + 2, 0);
                gaps._array[i + 1] = count;
                gaps._array[i + 2] = 0;
            }
        }
    }
    else  // gap bucket
    {
        gaps._array[i] += count;
    }

    // Adjust clipping information.
    gaps._clippingEndPos += count;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline typename Value<Gaps<TSequence, ArrayGaps> >::Type
value(Gaps<TSequence, ArrayGaps> const & gaps, TPosition clippedViewPos)
{
    if (isGap(gaps, clippedViewPos))
        return '-';
    else
        return value(source(gaps), toSourcePosition(gaps, clippedViewPos));
    return typename Value<Gaps<TSequence, ArrayGaps> >::Type();
}

// ----------------------------------------------------------------------------
// Function removeGaps()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition, typename TCount>
inline typename Size<Gaps<TSequence, ArrayGaps> >::Type
removeGaps(Gaps<TSequence, ArrayGaps> & gaps, TPosition clippedViewPos, TCount count)
{
    typedef Gaps<TSequence, ArrayGaps>     TGaps;
    typedef typename Position<TGaps>::Type TGapsPosition;
    typedef typename TGaps::TArray_        TArray;
    typedef typename TGaps::TArrayPos_     TArrayPos;
    typedef typename Value<TArray>::Type   TArrayValue;
    typedef typename Position<TSequence>::Type   TSeqPos;

    // Translate from clipped view position to unclipped view position.
    TGapsPosition pos = clippedViewPos + clippedBeginPosition(gaps);

    // Get index i of the according bucket and offset within bucket.
    SEQAN_ASSERT_GEQ(length(gaps._array), 2u);
    // Start at position 1 if there are no leading gaps.
    TArrayPos i = (gaps._array[0] == 0);
    TSeqPos offset = 0;
    for (TSeqPos counter = pos; counter > 0;)
    {
        SEQAN_ASSERT_LT(i, length(gaps._array));
        if (counter > gaps._array[i])
        {
            counter -= gaps._array[i];
            i += 1;
        }
        else
        {
            offset = counter;
            counter = 0;
        }
    }

    // Advance into next bucket if at end of current.
    if (offset > 0 && offset == gaps._array[i])
    {
        i += 1;
        offset = 0;
    }

    // If we are inside a non-gap bucket then we cannot remove any gaps.
    if (i % 2)
        return 0;

    // Otherwise, we can remove gaps right of the current position but not
    // more than there are.
    TSeqPos toRemove = count;
    if (toRemove > gaps._array[i] - offset)
        toRemove = gaps._array[i] - offset;
    gaps._array[i] -= toRemove;
    // In some cases, we remove the whole gap and merge the character buckets.
    if (gaps._array[i] == TArrayValue(0))
    {
        // No merging for leading and trailing gap.
        if (i != TArrayPos(0) && i != TArrayPos(length(gaps._array) - 1))
        {
            gaps._array[i - 1] += gaps._array[i + 1];
            erase(gaps._array, i, i + 2);
        }
    }

    // Also update the right clipping position.
    gaps._clippingEndPos -= toRemove;

    // Finally, return number of removed gaps.
    return toRemove;
}

// ----------------------------------------------------------------------------
// Function clearGaps()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline void
clearGaps(Gaps<TSequence, ArrayGaps> & gaps)
{
    _reinitArrayGaps(gaps);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline void
clear(Gaps<TSequence, ArrayGaps> & gaps)
{
    clear(gaps._source);
    clear(gaps._array);
    gaps._sourceBeginPos     = 0;
    gaps._sourceEndPos       = 0;
    gaps._clippingBeginPos   = 0;
    gaps._clippingEndPos     = 0;
    // cannot use clearGaps() here, since that calls value() on _source
    // which instates the Holder to Owner; we want it to be empty
}

// ----------------------------------------------------------------------------
// Function isGap()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline bool
isGap(Gaps<TSequence, ArrayGaps> const & gaps, TPosition clippedViewPos)
{
    typedef Gaps<TSequence, ArrayGaps>     TGaps;
    typedef typename Position<TGaps>::Type TGapsPosition;
    typedef typename TGaps::TArrayPos_     TArrayPos;
    typedef typename Position<TSequence>::Type TSeqPos;

    // Translate from clipped view position to unclipped view position.
    TGapsPosition pos = clippedViewPos + clippedBeginPosition(gaps);

    // Get index i of the according bucket and offset within bucket.
    SEQAN_ASSERT_GEQ(length(gaps._array), 2u);
    // Start at position 1 if there are no leading gaps.
    TArrayPos i = (gaps._array[0] == 0);
    TSeqPos offset = 0;
    for (TSeqPos counter = pos; counter > TSeqPos(0);)
    {
        SEQAN_ASSERT_LT(i, length(gaps._array));
        if (counter > gaps._array[i])
        {
            counter -= gaps._array[i];
            i += 1;
        }
        else
        {
            offset = counter;
            counter = 0;
        }
    }

    // Advance into next bucket if at end of current.
    if (offset > TSeqPos(0) && offset == gaps._array[i])
        i += 1;

    return !(i % 2);
}

// ----------------------------------------------------------------------------
// Function clearClipping()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline void
clearClipping(Gaps<TSequence, ArrayGaps> & gaps)
{
    typedef Gaps<TSequence, ArrayGaps>     TGaps;
    typedef typename TGaps::TArrayPos_     TArrayPos;

    gaps._sourceBeginPos = 0;
    gaps._sourceEndPos = length(value(gaps._source));
    gaps._clippingBeginPos = 0;
    gaps._clippingEndPos = 0;
    for (TArrayPos i = 0; i < length(gaps._array); ++i)
        gaps._clippingEndPos += gaps._array[i];

    SEQAN_ASSERT_LEQ(static_cast<int>(gaps._sourceEndPos), static_cast<int>(gaps._clippingEndPos));
}

// ----------------------------------------------------------------------------
// Function setClippedBeginPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline void
setClippedBeginPosition(Gaps<TSequence, ArrayGaps> & gaps, TPosition unclippedViewPosition)
{
    gaps._sourceBeginPos = toSourcePosition(gaps, unclippedViewPosition - clippedBeginPosition(gaps));
    gaps._clippingBeginPos = unclippedViewPosition;
}

// ----------------------------------------------------------------------------
// Function setClippedEndPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline void
setClippedEndPosition(Gaps<TSequence, ArrayGaps> & gaps, TPosition unclippedViewPosition)
{
    gaps._sourceEndPos = toSourcePosition(gaps, unclippedViewPosition - clippedBeginPosition(gaps));
    //if (isGap(gaps, unclippedViewPosition - clippedBeginPosition(gaps)))
    //    gaps._sourceEndPos += 1;
    gaps._clippingEndPos = unclippedViewPosition;
}

// ----------------------------------------------------------------------------
// Function clippedBeginPosition()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline typename Position<Gaps<TSequence, ArrayGaps> >::Type
clippedBeginPosition(Gaps<TSequence, ArrayGaps> const & gaps)
{
    return gaps._clippingBeginPos;
}

// ----------------------------------------------------------------------------
// Function clippedEndPosition()
// ----------------------------------------------------------------------------

template <typename TSequence>
inline typename Position<Gaps<TSequence, ArrayGaps> >::Type
clippedEndPosition(Gaps<TSequence, ArrayGaps> const & gaps)
{
    return gaps._clippingEndPos;
}

// ----------------------------------------------------------------------------
// Function setBeginPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline void
setBeginPosition(Gaps<TSequence, ArrayGaps> & gaps, TPosition sourcePosition)
{
    setClippedBeginPosition(gaps, toViewPosition(gaps, sourcePosition) + clippedBeginPosition(gaps));
}

// ----------------------------------------------------------------------------
// Function setEndPosition()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TPosition>
inline void
setEndPosition(Gaps<TSequence, ArrayGaps> & gaps, TPosition sourcePosition)
{
    setClippedEndPosition(gaps, toViewPosition(gaps, sourcePosition) + clippedBeginPosition(gaps));
}

// ----------------------------------------------------------------------------
// Function beginPosition()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We would rather like to have the const version only :(
template <typename TSequence>
inline typename Position<TSequence>::Type
beginPosition(Gaps<TSequence, ArrayGaps> const & gaps)
{
    return gaps._sourceBeginPos;
}

template <typename TSequence>
inline typename Position<TSequence>::Type
beginPosition(Gaps<TSequence, ArrayGaps> & gaps)
{
    return gaps._sourceBeginPos;
}

// ----------------------------------------------------------------------------
// Function endPosition()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We would rather like to have the const version only :(
template <typename TSequence>
inline typename Position<TSequence>::Type
endPosition(Gaps<TSequence, ArrayGaps> const & gaps)
{
    return gaps._sourceEndPos;
}

template <typename TSequence>
inline typename Position<ArrayGaps>::Type
endPosition(Gaps<TSequence, ArrayGaps> & gaps)
{
    return gaps._sourceEndPos;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GAPS_ARRAY_H_
