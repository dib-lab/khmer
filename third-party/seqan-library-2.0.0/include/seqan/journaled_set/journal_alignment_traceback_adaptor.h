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
// Adapts traceback to a journaled string using SortedArray specialization.
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNAL_ALIGNMENT_TRACEBACK_ADAPTOR_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNAL_ALIGNMENT_TRACEBACK_ADAPTOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _adaptTraceSegmentsTo()           [String<Journaled<SortedArray> >]
// ----------------------------------------------------------------------------


template <typename TValue, typename THostSpec, typename TBuffSpec, typename TReference,
          typename TSource, typename TPosition, typename TSize, typename TStringSpec>
inline void
_adaptTraceSegmentsTo(String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > & targetJournal,
                      TReference const & /*reference*/, //TODO(rmaerker): What if the TSequenceH is a journal string.
                      TSource const & source,
                      String<TraceSegment_<TPosition, TSize>, TStringSpec> const & traceSegments)
{

    typedef String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > TJournalString;
    typedef typename JournalType<TJournalString>::Type TJournalEntries;
    typedef typename Value<TJournalEntries>::Type TJournalEntry;

//    typedef typename Position<TReference>::Type TPhysicalPosition;
    typedef typename Position<TSource>::Type TVirtualPosition;

    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef typename Iterator<String<TTraceSegment, TStringSpec> const>::Type TTraceIterator;

    SEQAN_ASSERT_NOT_MSG(empty(host(targetJournal)), "No reference sequence set!");

    // TODO(rmaerker): Hot-fix for new alignment module.
    clear(targetJournal._journalEntries._journalNodes);
    clear(targetJournal._insertionBuffer);
    targetJournal._journalEntries._originalStringLength = length(host(targetJournal));

    TVirtualPosition virtualPos = 0;

    TTraceIterator traceSegIter = end(traceSegments, Standard()) - 1;
    TTraceIterator traceSegIterEnd = begin(traceSegments, Standard()) - 1;

    for (; traceSegIter != traceSegIterEnd; --traceSegIter)
    {
        switch (value(traceSegIter)._traceValue)
        {
        case TraceBitMap_::DIAGONAL: //matching area
        {
            appendValue(targetJournal._journalEntries._journalNodes, TJournalEntry(SOURCE_ORIGINAL, value(traceSegIter)._horizontalBeginPos,
                                                                                   virtualPos, value(traceSegIter)._horizontalBeginPos, value(traceSegIter)._length));
            virtualPos += value(traceSegIter)._length;
            break;
        }

        case TraceBitMap_::VERTICAL: //insertion
        {
            appendValue(targetJournal._journalEntries._journalNodes, TJournalEntry(SOURCE_PATCH, length(targetJournal._insertionBuffer),
                                                                                   virtualPos, 0, value(traceSegIter)._length));
            append(targetJournal._insertionBuffer,
                   infix(source, value(traceSegIter)._verticalBeginPos, value(traceSegIter)._verticalBeginPos +
                         value(traceSegIter)._length));
            virtualPos += value(traceSegIter)._length;
            break;
        }

        default:
            break;  // Otherwise we are in a deletion and do nothing.
        }
    }
    _setLength(targetJournal, virtualPos);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNAL_ALIGNMENT_TRACEBACK_ADAPTOR_H_
