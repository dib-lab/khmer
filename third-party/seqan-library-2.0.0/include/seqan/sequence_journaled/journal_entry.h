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
// Code for journal entries.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRY_H_
#define SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRY_H_

namespace seqan {

// ============================================================================
// Enums, Classes
// ============================================================================

enum SegmentSource
{
    SOURCE_NULL,
    SOURCE_ORIGINAL,
    SOURCE_PATCH
};


template <typename TPos_, typename TSize_>
struct JournalEntry
{
    typedef TPos_ TPos;
    typedef TSize_ TSize;

    // Flag for where the segment comes from.
    SegmentSource segmentSource;
    // Position in the original string or the insertion buffer,
    // depending on segmentSource.
    TPos physicalPosition;
    // Position in the virtual string.
    TPos virtualPosition;
    // Physical position in the host string.  For SOURCE_PATCH entries, this is the position of the closest
    // SOURCE_ORIGINAL entry left of it.  If there is no such entry then this is 0.  Unused and always 0 for unbalanced
    // tree.
    TPos physicalOriginPosition;
    // Length of the segment.
    TSize length;

    JournalEntry()
            : segmentSource(SOURCE_NULL),
              physicalPosition(0),
              virtualPosition(0),
              physicalOriginPosition(0),
              length(0)
    {}

    JournalEntry(SegmentSource const & segmentSource,
                 TPos physicalPosition,
                 TPos virtualPosition,
                 TPos physicalHostPosition,
                 TSize length)
            : segmentSource(segmentSource),
              physicalPosition(physicalPosition),
              virtualPosition(virtualPosition),
              physicalOriginPosition(physicalHostPosition),
              length(length)
    {}
};


template <typename TPos, typename TSize>
struct JournalEntryLtByVirtualPos
{
    bool operator()(JournalEntry<TPos, TSize> const & a,
                    JournalEntry<TPos, TSize> const & b) const
    {
        return a.virtualPosition < b.virtualPosition;
    }
};


// TODO (rmaerker): Note that this comparison only makes sense if the nodes are sorted
// in ascending order according to their physical position - this might be violated
// if rearrangements are included.
template <typename TPos, typename TSize>
struct JournalEntryLtByPhysicalOriginPos
{
    bool operator()(JournalEntry<TPos, TSize> const & a, JournalEntry<TPos, TSize> const & b) const
    {
        return a.physicalOriginPosition < b.physicalOriginPosition;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TPos, typename TSize>
struct Size<JournalEntry<TPos, TSize> >
{
    typedef TSize Type;
};

template <typename TPos, typename TSize>
struct Size<JournalEntry<TPos, TSize> const>
        : Size<JournalEntry<TPos, TSize> > {};

template <typename TPos, typename TSize>
struct Position<JournalEntry<TPos, TSize> >
{
    typedef TPos Type;
};

template <typename TPos, typename TSize>
struct Position<JournalEntry<TPos, TSize> const>
        : Position<JournalEntry<TPos, TSize> > {};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TPos, typename TSize>
TStream & operator<<(TStream & stream, JournalEntry<TPos, TSize> const & entry)
{
    return stream << "{segmentSource=" << entry.segmentSource
                  << ", virtualPosition=" << entry.virtualPosition
                  << ", physicalPosition=" << entry.physicalPosition
                  << ", physicalOriginPosition=" << entry.physicalOriginPosition
                  << ", length=" << entry.length
                  << "}";
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRY_H_

