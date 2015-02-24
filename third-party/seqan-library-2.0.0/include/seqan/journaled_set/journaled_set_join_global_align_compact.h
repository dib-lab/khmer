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

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_GLOBAL_ALIGN_COMPACT_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_GLOBAL_ALIGN_COMPACT_H_

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
// Function _joinInternal()                     [GlobalAlign<JournaledCompact>]
// ----------------------------------------------------------------------------

// TODO(rmaerker): Change the orientation of reference and journal.
template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TJournalString2>
inline void
_joinInternal(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journal,
              StringSet<TJournalString2, Owner<JournaledSet> > const & journalSet,
              JoinConfig<GlobalAlign<JournaledCompact> > const & joinConfig)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;

    // TODO(rmaerker): Check the correct behavior here.
    TJournalString tmpJournal;
    setHost(tmpJournal, value(journalSet._globalRefHolder));

    if (isBandSet(joinConfig))
        globalAlignment(tmpJournal, host(journalSet), journal, scoringScheme(joinConfig), joinConfig._alignConfig,
                        lowerDiagonal(joinConfig), upperDiagonal(joinConfig));
    else
        globalAlignment(tmpJournal, host(journalSet), journal, scoringScheme(joinConfig), joinConfig._alignConfig);

    // Apply the alignment to the journal.
    set(journal, tmpJournal);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_GLOBAL_ALIGN_COMPACT_H_
