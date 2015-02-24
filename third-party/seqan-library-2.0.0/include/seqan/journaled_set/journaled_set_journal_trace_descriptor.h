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

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOURNAL_TRACE_DESCRIPTOR_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOURNAL_TRACE_DESCRIPTOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Enum TraceDirection
// ----------------------------------------------------------------------------

enum TraceDirection
{
        DIAGONAL,
        HORIZONTAL,
        VERTICAL
};

// ----------------------------------------------------------------------------
// Class JournalTraceBuffer
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
class JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > >
{
public:
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString;
    typedef typename Position<TString>::Type TPos;
    typedef typename Size<TString>::Type TSize;
    typedef JournalEntry<TPos, TSize> TJournalNode;

    // Stores the journal operations received from the trace back in reverse order.
    String<TJournalNode> revSortedOperation_;
    String<TValue, TBuffSpec> insertionBuffer_;

    JournalTraceBuffer()
    {
    }

    template <typename TPos>
    inline typename Reference<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >::Type
    operator[](TPos const pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>::Type
    operator[](TPos const pos) const
    {
        return value(*this, pos);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Value<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
{
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
    typedef typename Position<TString_>::Type TPos_;
    typedef typename Size<TString_>::Type TSize_;
    typedef JournalEntry<TPos_, TSize_> Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Value<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
{
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
    typedef typename Position<TString_>::Type TPos_;
    typedef typename Size<TString_>::Type TSize_;
    typedef JournalEntry<TPos_, TSize_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Reference<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
{
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
    typedef typename Value<JournalTraceBuffer<TString_> >::Type & Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Reference<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
{
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
    typedef typename Value<JournalTraceBuffer<TString_> >::Type const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Size<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
{
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
    typedef typename Size<TString_>::Type Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Size<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
    : Size<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Position<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
{
    typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
    typedef typename Position<TString_>::Type Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Position<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
    : Position<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > > {};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Iterator<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > >, Standard >
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString_;
    typedef typename TJournalString_::TJournalEntry TJournalEntry_;
    typedef typename Iterator< String<TJournalEntry_>, Standard >::Type Type;
};

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
struct Iterator<JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const, Standard >
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString_;
    typedef typename TJournalString_::TJournalEntry TJournalEntry_;
    typedef typename Iterator< String<TJournalEntry_> const, Standard >::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Iterator<JournalTraceBuffer<TString>, Standard >::Type
begin(JournalTraceBuffer<TString> & me)
{
    SEQAN_CHECKPOINT;

    return begin(getTrace(me));
}

template <typename TString>
inline typename Iterator<JournalTraceBuffer<TString> const, Standard >::Type
begin(JournalTraceBuffer<TString> const & me)
{
    return begin(getTrace(me));
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Iterator<JournalTraceBuffer<TString>, Standard >::Type
end(JournalTraceBuffer<TString> & me)
{
    return end(getTrace(me));
}

template <typename TString>
inline typename Iterator<JournalTraceBuffer<TString> const, Standard >::Type
end(JournalTraceBuffer<TString> const & me)
{
    return end(getTrace(me));
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TJournal>
TStream &
operator <<(TStream & stream,
            JournalTraceBuffer<TJournal> const & obj)
{
    for (unsigned int i = 0; i < length(obj);++i)
    {
        if (obj[i].segmentSource == SOURCE_PATCH)
        {
            stream << obj[i] << " " <<infix(obj.insertionBuffer_,obj[i].physicalPosition, obj[i].physicalPosition + obj[i].length)<< std::endl;
        } else
        {
            stream << obj[i] << std::endl;
        }
    }
    return stream;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TJournal>
inline typename Size<JournalTraceBuffer<TJournal> >::Type
length(JournalTraceBuffer<TJournal> const & me)
{
    return length(me.revSortedOperation_);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TJournalString>
inline void
clear(JournalTraceBuffer<TJournalString> & traceBuff)
{
    clear(traceBuff.revSortedOperation_);
    clear(traceBuff.insertionBuffer_);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------
template <typename TJournalString>
inline typename Reference<JournalTraceBuffer<TJournalString> >::Type
value(JournalTraceBuffer<TJournalString> & me,
      typename Position<JournalTraceBuffer<TJournalString> >::Type const pos)
{
    SEQAN_CHECKPOINT;
    return me.revSortedOperation_[pos];
}

template <typename TJournalString>
inline typename Reference<JournalTraceBuffer<TJournalString> const>::Type
value(JournalTraceBuffer<TJournalString> const & me,
      typename Position<JournalTraceBuffer<TJournalString> >::Type const pos)
{
    return me.revSortedOperation_[pos];
}

// ----------------------------------------------------------------------------
// Function getTrace()
// ----------------------------------------------------------------------------

template <typename TString>
inline String <typename Value<JournalTraceBuffer<TString> >::Type> &
getTrace(JournalTraceBuffer<TString > & me)
{
    return me.revSortedOperation_;
}


template <typename TString>
inline String<typename Value<JournalTraceBuffer<TString> >::Type> const &
getTrace(JournalTraceBuffer<TString> const & me)
{
    return me.revSortedOperation_;
}

// ----------------------------------------------------------------------------
// Function reverseTrace()
// ----------------------------------------------------------------------------

// TODO(rmaerker): Use modifier instead.
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
void inline
reverseTrace(JournalTraceBuffer<String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me)
{
    reverse(me.revSortedOperation_);
}

// ----------------------------------------------------------------------------
// Function getTraceReverse()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
inline String <typename Value<JournalTraceBuffer<String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >::Type>
getTraceReverse(JournalTraceBuffer<String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me)
{
    SEQAN_CHECKPOINT;

    typedef String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournal;
    typedef typename Value<JournalTraceBuffer<TJournal> >::Type  TEntry;
    typedef String<TEntry> TJournalEntries;

    TJournalEntries cpy(me.revSortedOperation_);
    reverse(cpy);
    return cpy;
}

template <typename TString>
inline String <typename Value<JournalTraceBuffer<TString> const>::Type> const
getTraceReverse(JournalTraceBuffer<TString> const & me)
{
    SEQAN_CHECKPOINT;
    return getTraceReverse(const_cast<JournalTraceBuffer<TString> &>(me));
}

// ----------------------------------------------------------------------------
// Function getInsertionBuffer()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
inline String<TValue, TBuffSpec> &
getInsertionBuffer(JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me)
{
    SEQAN_CHECKPOINT;
    return me.insertionBuffer_;
}

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
inline String<TValue, TBuffSpec> const &
getInsertionBuffer(JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const & me)

{
    SEQAN_CHECKPOINT;
    return me.insertionBuffer_;
}


// ----------------------------------------------------------------------------
// Function _applyTraceOperations()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
void
_applyTraceOperations(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journal,
                      String<TValue, THostSpec> const & newHost,
                      JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const & traceDescr)
{
    typedef String<TValue, THostSpec> THost;

    setValue(journal._holder, const_cast<THost &>(newHost));
    reinit(journal._journalEntries, length(newHost));
    _constructAndSetJournalTree(journal, getTrace(traceDescr), getInsertionBuffer(traceDescr));
}

// ----------------------------------------------------------------------------
// Function _alignTracePrint()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TSequenceH,
    typename TSequenceV, typename TId, typename TPos, typename TTraceValue>
inline void
_alignTracePrint(JournalTraceBuffer<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me,
         TSequenceH const & /*seqH*/,
         TSequenceV const & seqV,
         TId,
         TPos const pos1,
         TId,
         TPos const pos2,
         TPos const segLen,
         TTraceValue const tv)
 {
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TString;
    typedef typename Value<JournalTraceBuffer<TString> >::Type TEntryString;
    typedef typename Value<TEntryString>::Type TJournalEntry;

    enum SegmentSource segmentSrc = SOURCE_NULL;
    TPos physicalPos = 0;
    TPos physicalOriginPos = 0;

    if (segLen == 0)
    {
        return;
    }
    switch (tv)
    {
    case DIAGONAL://matching area
    {
        segmentSrc = SOURCE_ORIGINAL;
        physicalPos = pos1;
        physicalOriginPos = pos1;
        break;
    }
    case VERTICAL://insertion
    {
        segmentSrc = SOURCE_PATCH;
        physicalPos = length(me.insertionBuffer_);
        append(me.insertionBuffer_, infix(seqV,pos2, pos2 + segLen));
        break;
    }
    case HORIZONTAL://deletion - nothing to be done here
        return;
        break;
    }
    // TODO(rmaerker): Change computation of physicalOrigin position when porting to new alignment module.
    appendValue(getTrace(me), TJournalEntry(segmentSrc, physicalPos, pos2, physicalOriginPos, segLen));
 }

// ----------------------------------------------------------------------------
// Function _constructAndSetJournalTree()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TBuffSpec, typename TCargo, typename THostSpec2>
void
_constructAndSetJournalTree(String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > & journalSeq,
                            String<TCargo, THostSpec2> const & cargoArray,
                            String<TValue, TBuffSpec> const & insertionBuffer)
{
    assign(journalSeq._journalEntries._journalNodes, cargoArray);
    assign(journalSeq._insertionBuffer, insertionBuffer);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOURNAL_TRACE_DESCRIPTOR_H_
