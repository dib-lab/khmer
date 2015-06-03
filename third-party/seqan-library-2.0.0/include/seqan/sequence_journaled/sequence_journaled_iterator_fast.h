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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements a fast version of the iterator over journal strings.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEQUENCE_JOURNALED_SEQUENCE_JOURNALED_ITERATOR_FAST_H_
#define INCLUDE_SEQAN_SEQUENCE_JOURNALED_SEQUENCE_JOURNALED_ITERATOR_FAST_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TJournaledString>
class Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> >
{
public:
    typedef Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > TIterator;
    typedef typename JournalType<TJournaledString>::Type TJournalEntries;
    // We need a rooted iterator for iterating the journal tree since we need atEnd().
    typedef typename Iterator<TJournalEntries, Standard>::Type TJournalEntriesIterator;
    typedef typename Host<TJournaledString>::Type THost;
    typedef typename Iterator<THost, Standard>::Type TSegmentIterator;

    // The journal string we iterate over.
    TJournaledString * _journalStringPtr;
    // Iterator over the segments in the journal tree.
    TJournalEntriesIterator _journalEntriesIterator;
    // Begin and end iterator in the host string of the journal string.
    TSegmentIterator _segmentBegin;
    TSegmentIterator _segmentEnd;
    // Current iterator in the host segment.
    TSegmentIterator _currentSegmentIt;

    Iter() :
        _journalStringPtr(), _journalEntriesIterator(), _segmentBegin(), _segmentEnd(),
        _currentSegmentIt()
    {}

    Iter(TIterator const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _segmentBegin(other._segmentBegin),
              _segmentEnd(other._segmentEnd),
              _currentSegmentIt(other._currentSegmentIt)
    {
        SEQAN_CHECKPOINT;
    }

    Iter(typename IterComplementConst<TIterator>::Type const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _segmentBegin(other._segmentBegin),
              _segmentEnd(other._segmentEnd),
              _currentSegmentIt(other._currentSegmentIt)
    {
        SEQAN_CHECKPOINT;
    }

    Iter & operator=(TIterator const & other)
    {
        if (this != &other)
        {
            _journalStringPtr = other._journalStringPtr;
            _journalEntriesIterator = other._journalEntriesIterator;
            _segmentBegin = other._segmentBegin;
            _segmentEnd = other._segmentEnd;
            _currentSegmentIt = other._currentSegmentIt;
        }
        return *this;
    }

    Iter & operator=(typename IterComplementConst<TIterator>::Type const & other)
    {
        if (this != &other)
        {
            _journalStringPtr = other._journalStringPtr;
            _journalEntriesIterator = other._journalEntriesIterator;
            _segmentBegin = other._segmentBegin;
            _segmentEnd = other._segmentEnd;
            _currentSegmentIt = other._currentSegmentIt;
        }
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
// Function _initJournaledStringIteratorEnd()
// ----------------------------------------------------------------------------

template <typename TJournaledString>
inline
void
_initJournaledStringIteratorEnd(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator)
{
    iterator._journalEntriesIterator = end(iterator._journalStringPtr->_journalEntries, Standard()) - 1;
    _updateSegmentIteratorsLeft(iterator);
    ++iterator._currentSegmentIt;
}

// ----------------------------------------------------------------------------
// Function _updateSegmentIterators()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Rename to _updateSegmentIteratorsRight().
template <typename TJournaledString>
inline
void
_updateSegmentIterators(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator)
{
    if (atEnd(iterator._journalEntriesIterator, _journalEntries(container(iterator))))
    {
        _initJournaledStringIteratorEnd(iterator);
        return;
    }

    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL)
        iterator._segmentBegin = begin(host(*iterator._journalStringPtr), Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
    else
        iterator._segmentBegin = begin(iterator._journalStringPtr->_insertionBuffer, Standard()) + value(iterator._journalEntriesIterator).physicalPosition;

    iterator._segmentEnd = iterator._segmentBegin + value(iterator._journalEntriesIterator).length;
    iterator._currentSegmentIt = iterator._segmentBegin;
}

// ----------------------------------------------------------------------------
// Function _updateSegmentIteratorsLeft()
// ----------------------------------------------------------------------------

template <typename TJournaledString>
inline
void
_updateSegmentIteratorsLeft(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator)
{
    if (iterator._journalEntriesIterator + 1 == begin(_journalEntries(container(iterator)), Standard()))
    {
        goBegin(iterator);
        return;
    }

    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL)
        iterator._segmentBegin = begin(host(*iterator._journalStringPtr), Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
    else
        iterator._segmentBegin = begin(iterator._journalStringPtr->_insertionBuffer, Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
    iterator._segmentEnd = iterator._segmentBegin + value(iterator._journalEntriesIterator).length;
    iterator._currentSegmentIt = iterator._segmentEnd - 1;
}

// ----------------------------------------------------------------------------
// Function _localEntryPosition
// ----------------------------------------------------------------------------

// TODO(rmaerker): Write documentation!
template <typename TJournaledString>
inline typename Position<Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const>::Type
_localEntryPosition(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const & iterator)
{
    return iterator._currentSegmentIt - iterator._segmentBegin;
}

// ----------------------------------------------------------------------------
// Function setPosition
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TPosition>
inline void
setPosition(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & me,
            TPosition pos)
{
    SEQAN_ASSERT_GEQ(pos, static_cast<TPosition>(0));

    // Handle case where pos points behind the container.
    if (pos >= static_cast<TPosition>(length(container(me))))
        goEnd(me);

    // Use binary search to find corresponding node.
    me._journalEntriesIterator = findInJournalEntries(container(me)._journalEntries, pos);

    int offset = pos - value(me._journalEntriesIterator).virtualPosition;  // The offset within the current node.
    _updateSegmentIterators(me);
    me._currentSegmentIt += offset;
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

// getValue
template <typename TJournaledString>
inline
typename GetValue<Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > >::Type
getValue(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const & iterator)
{
    return getValue(iterator._currentSegmentIt);
}

template <typename TJournaledString>
inline
typename GetValue<Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > >::Type
getValue(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator)
{
    return getValue(iterator._currentSegmentIt);
}

// ----------------------------------------------------------------------------
// Function operator++()
// ----------------------------------------------------------------------------

template <typename TJournaledString>
inline
Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > &
operator++(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator)
{
    ++iterator._currentSegmentIt;
    if (iterator._currentSegmentIt == iterator._segmentEnd)
    {
        ++iterator._journalEntriesIterator;
        _updateSegmentIterators(iterator);
    }
    return iterator;
}

// ----------------------------------------------------------------------------
// Function opertor+=()
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TLen>
inline
Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > &
operator+=(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator,
           TLen len_)
{
    SEQAN_ASSERT_GEQ(len_, static_cast<TLen>(0));

    TLen remaining = iterator._segmentEnd - iterator._currentSegmentIt;
    while (len_ > 0  && remaining != 0)
    {
        SEQAN_ASSERT_GT(remaining, static_cast<TLen>(0));
        if (len_ >= remaining)
        {
            len_ -= remaining;
            ++iterator._journalEntriesIterator;
            _updateSegmentIterators(iterator);
            remaining = iterator._segmentEnd - iterator._currentSegmentIt;
        }
        else
        {
            iterator._currentSegmentIt += len_;
            len_ = 0;
        }
    }
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator--()                                               [prefix]
// ----------------------------------------------------------------------------

template <typename TJournaledString>
inline
Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > &
operator--(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator)
{
    if (iterator._currentSegmentIt == iterator._segmentBegin)
    {
        --iterator._journalEntriesIterator;
        _updateSegmentIteratorsLeft(iterator);
        return iterator;
    }

    --iterator._currentSegmentIt;
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TJournaledString, typename TLen>
inline
Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > &
operator-=(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > & iterator,
            TLen len_)
{
    typedef Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > TIterator;
    typedef typename Size<TIterator>::Type TSize;

    SEQAN_ASSERT_GEQ(len_, static_cast<TLen>(0));
    TSize len = len_;

    while (len > 0 )
    {
        TSize relNodePos = _localEntryPosition(iterator);
        if (len > relNodePos)
        {
            len -= (relNodePos + 1);
            --iterator._journalEntriesIterator;
            _updateSegmentIteratorsLeft(iterator);
        }
        else
        {
            iterator._currentSegmentIt -= len;
            len = 0;
        }
    }
    return iterator;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TJournaledString>
inline
typename Difference<TJournaledString>::Type
operator-(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const & it1,
          Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const & it2)
{
    // First, handle the cases where it1 or it2 are at the end.
    bool it1AtEnd = atEnd(it1._journalEntriesIterator, _journalEntries(container(it1)));
    bool it2AtEnd = atEnd(it2._journalEntriesIterator, _journalEntries(container(it2)));
    if (it1AtEnd && it2AtEnd)
    {
        return 0;
    }
    else if (it1AtEnd)
    {
        SEQAN_ASSERT_LT(value(it2._journalEntriesIterator).virtualPosition + _localEntryPosition(it2), length(*it1._journalStringPtr));
        return length(*it1._journalStringPtr) - (value(it2._journalEntriesIterator).virtualPosition + _localEntryPosition(it2));
    }
    else if (it2AtEnd)
    {
        SEQAN_ASSERT_LT((value(it1._journalEntriesIterator).virtualPosition + _localEntryPosition(it1)), length(*it1._journalStringPtr));
        return (value(it1._journalEntriesIterator).virtualPosition + _localEntryPosition(it1)) - length(*it1._journalStringPtr);
    }
    // Otherwise, we can simply subtract the virtual positions.
    return (value(it1._journalEntriesIterator).virtualPosition + _localEntryPosition(it1)) - (value(it2._journalEntriesIterator).virtualPosition + _localEntryPosition(it2));
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TJournaledString>
inline
bool
operator==(Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const & a,
           Iter<TJournaledString, JournaledStringIterSpec<CommonSegmentIterator> > const & b)
{
    SEQAN_CHECKPOINT;
    if (a._journalEntriesIterator != b._journalEntriesIterator)
        return false;
    if (a._currentSegmentIt != b._currentSegmentIt)
        return false;
    return true;
}

}

#endif // INCLUDE_SEQAN_SEQUENCE_JOURNALED_SEQUENCE_JOURNALED_ITERATOR_FAST_H_
