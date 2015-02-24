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

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_GLOBAL_ALIGN_MANHATTEN_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_GLOBAL_ALIGN_MANHATTEN_H_

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
// Function _joinInternal()                   [GlobalAlign<JournaledManhatten>]
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TJournalString2>
inline void
_joinInternal(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journal,
              StringSet<TJournalString2, Owner<JournaledSet> > const & journalSet,
              JoinConfig<GlobalAlign<JournaledManhatten> > const &)
{

    typedef String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
    typedef typename TJournalString::TJournalEntry TJounralEntry;
//    typedef typename Size<TJournalString>::Type TSize;

    JournalTraceBuffer<TJournalString> traceDescriptor;
    std::stringstream stream;
    stream << journal;
    appendValue(getTrace(traceDescriptor), TJounralEntry(SOURCE_PATCH, 0, 0, 0,length(journal)));
    append(getInsertionBuffer(traceDescriptor), stream.str());
    _applyTraceOperations(journal, host(journalSet), traceDescriptor);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_GLOBAL_ALIGN_MANHATTEN_H_
