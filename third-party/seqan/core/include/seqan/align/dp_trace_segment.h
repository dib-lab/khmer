// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
//
// The TraceSegment structure is used to store the traceback in a common
// structure such that we can easiely adapt them afterwards in the
// user-defined structure, such as Align or AlignmentGraph objects.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_TRACE_SEGMENT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_TRACE_SEGMENT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class TraceSegment
// ----------------------------------------------------------------------------

// TraceSegments are used as a common interface to all structures that can represent an alignment.
//
// See alignment_dp_traceback_adaptor.h to find methods for adaption.
template <typename TPosition, typename TSize>
class TraceSegment_
{
public:
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    TPosition _horizontalBeginPos;      // the begin position in horizontal dimension
    TPosition _verticalBeginPos;        // the begin position in vertical dimension
    TSize _length;                      // the length of the segment
    TTraceValue _traceValue;            // the trace direction

    TraceSegment_() :
        _horizontalBeginPos(0), _verticalBeginPos(0), _length(0), _traceValue(+TraceBitMap_::NONE){}

    TraceSegment_(TraceSegment_ const & other) :
        _horizontalBeginPos(other._horizontalBeginPos),
        _verticalBeginPos(other._verticalBeginPos),
        _length(other._length),
        _traceValue(other._traceValue) {}

    TraceSegment_(TPosition const & horizontalBeginPos, TPosition const & verticalBeginPos, TSize const & length,
                  TTraceValue const & traceValue) :
        _horizontalBeginPos(horizontalBeginPos),
        _verticalBeginPos(verticalBeginPos),
        _length(length),
        _traceValue(traceValue) {}

    TraceSegment_ &
    operator=(TraceSegment_ const & other)
    {
        if (this != &other)
        {
            _horizontalBeginPos = other._horizontalBeginPos;
            _verticalBeginPos = other._verticalBeginPos;
            _length = other._length;
            _traceValue = other._traceValue;
        }
        return *this;
    }

};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
struct Position<TraceSegment_<TPosition, TSize> >
{
    typedef TPosition Type;
};

template <typename TPosition, typename TSize>
struct Position<TraceSegment_<TPosition, TSize> const>:
    Position<TraceSegment_<TPosition, TSize> >{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
struct Size<TraceSegment_<TPosition, TSize> >
{
    typedef TSize Type;
};

template <typename TPosition, typename TSize>
struct Size<TraceSegment_<TPosition, TSize> const>:
    Size<TraceSegment_<TPosition, TSize> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _getBeginHorizontal()
// ----------------------------------------------------------------------------

// The begin position of the segment in horizontal dimension.
template <typename TPosition, typename TSize>
inline TPosition
_getBeginHorizontal(TraceSegment_<TPosition, TSize> const & traceSegment)
{
    return traceSegment._horizontalBeginPos;
}

// ----------------------------------------------------------------------------
// Function _getBeginVertical()
// ----------------------------------------------------------------------------

// The begin position of the segment in vertical dimension.
template <typename TPosition, typename TSize>
inline TPosition
_getBeginVertical(TraceSegment_<TPosition, TSize> const & traceSegment)
{
    return traceSegment._verticalBeginPos;
}

// ----------------------------------------------------------------------------
// Function _getEndHorizontal()
// ----------------------------------------------------------------------------

// The end position of the segment in horizontal dimension.
template <typename TPosition, typename TSize>
inline TPosition
_getEndHorizontal(TraceSegment_<TPosition, TSize> const & traceSegment)
{
    if (traceSegment._traceValue & (TraceBitMap_::HORIZONTAL | TraceBitMap_::DIAGONAL))
    {
        return traceSegment._horizontalBeginPos + traceSegment._length;
    }
    return traceSegment._horizontalBeginPos;
}

// ----------------------------------------------------------------------------
// Function _getEndVertical()
// ----------------------------------------------------------------------------

// The end position of the segment in vertical dimension.
template <typename TPosition, typename TSize>
inline TPosition
_getEndVertical(TraceSegment_<TPosition, TSize> const & traceSegment)
{
    if (traceSegment._traceValue & (TraceBitMap_::VERTICAL | TraceBitMap_::DIAGONAL))
    {
        return traceSegment._verticalBeginPos + traceSegment._length;
    }
    return traceSegment._verticalBeginPos;
}

// ----------------------------------------------------------------------------
// Function _getTraceValue()
// ----------------------------------------------------------------------------

// The end position of the segment in vertical dimension.
template <typename TPosition, typename TSize>
inline typename TraceSegment_<TPosition, TSize>::TTraceValue
_getTraceValue(TraceSegment_<TPosition, TSize> const & traceSegment)
{
    return traceSegment._traceValue;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// The length of the segment.
template <typename TPosition, typename TSize>
inline TSize
length(TraceSegment_<TPosition, TSize> const & traceSegment)
{
    return traceSegment._length;
}

// ----------------------------------------------------------------------------
// Function _setLength()
// ----------------------------------------------------------------------------

// The length of the segment.
template <typename TPosition, typename TSize>
inline void
_setLength(TraceSegment_<TPosition, TSize> & traceSegment, TSize newLength)
{
    traceSegment._length = newLength;
}



// ----------------------------------------------------------------------------
// Function _translateTraceValue()
// ----------------------------------------------------------------------------

// Translates the trace value into a human-readable format.
//
// Note, used for debugging reasons only.
template <typename TTraceValue>
String<char> _translateTraceValue(TTraceValue const & traceValue)
{
    String<char> transcript;

    if ((traceValue & TraceBitMap_::DIAGONAL) == TraceBitMap_::DIAGONAL)
        append(transcript, 'D');
    if ((traceValue & TraceBitMap_::VERTICAL) == TraceBitMap_::VERTICAL)
        append(transcript, 'V');
    if ((traceValue & TraceBitMap_::HORIZONTAL) == TraceBitMap_::HORIZONTAL)
        append(transcript, 'H');
    if ((traceValue & TraceBitMap_::VERTICAL_OPEN) == TraceBitMap_::VERTICAL_OPEN)
        append(transcript, 'v');
    if ((traceValue & TraceBitMap_::HORIZONTAL_OPEN) == TraceBitMap_::HORIZONTAL_OPEN)
        append(transcript, 'h');
    if ((traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX) == TraceBitMap_::MAX_FROM_VERTICAL_MATRIX)
        append(transcript, '|');
    if ((traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX) == TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX)
        append(transcript, '-');

    if ((traceValue) == TraceBitMap_::NONE)
        append(transcript, '0');
    return transcript;
}

// ----------------------------------------------------------------------------
// Function oerpator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TSize, typename TPosition>
TStream & operator<<(TStream & stream, TraceSegment_<TSize, TPosition> const & traceSegment)
{
    stream << _translateTraceValue(traceSegment._traceValue) << "-";
    stream << "(" << traceSegment._horizontalBeginPos << ", " << traceSegment._verticalBeginPos << ", " <<
    traceSegment._length << ")";
    return stream;
}

// ----------------------------------------------------------------------------
// Function oerpator==()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline bool operator==(TraceSegment_<TPosition, TSize> const & left, TraceSegment_<TPosition, TSize> const & right)
{
    if (left._horizontalBeginPos != right._horizontalBeginPos)
        return false;

    if (left._verticalBeginPos != right._verticalBeginPos)
        return false;

    if (left._length != right._length)
        return false;

    if (left._traceValue != right._traceValue)
        return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function oerpator!=()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline bool operator!=(TraceSegment_<TPosition, TSize> const & left, TraceSegment_<TPosition, TSize> const & right)
{
    return !(left == right);
}

// ----------------------------------------------------------------------------
// Function recordSegment()
// ----------------------------------------------------------------------------

// Records a segment given the horizontal and vertical begin position, the length,
// and the corrsponding trace value.
//
// The first parameter is the container the segment is recorded to.
template <typename TTraceSegments, typename TPositionH, typename TPositionV, typename TSize, typename TTraceValue>
inline void _recordSegment(TTraceSegments & traceSegments,
                           TPositionH const & horizontalBeginPos,
                           TPositionV const & verticalBeginPos,
                           TSize const & segmentLength,
                           TTraceValue const & traceValue)
{
    typedef typename Value<TTraceSegments>::Type TTraceSegment;

    if (segmentLength == 0)
        return;  // we don't store empty segments

    if (traceValue & TraceBitMap_::DIAGONAL)
        appendValue(traceSegments, TTraceSegment(horizontalBeginPos, verticalBeginPos, segmentLength, +TraceBitMap_::DIAGONAL));
    else if (traceValue & TraceBitMap_::VERTICAL)
        appendValue(traceSegments, TTraceSegment(horizontalBeginPos, verticalBeginPos, segmentLength, +TraceBitMap_::VERTICAL));
    else if (traceValue & TraceBitMap_::HORIZONTAL)
        appendValue(traceSegments, TTraceSegment(horizontalBeginPos, verticalBeginPos, segmentLength, +TraceBitMap_::HORIZONTAL));
    // everything else is not tracked.
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_TRACE_SEGMENT_H_
