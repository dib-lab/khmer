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

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_H_

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

/*!
 * @fn JournaledSet#join
 * @brief Joins a @link JournaledString @endlink to a @link JournaledSet @endlink by computing and journaling
 *        differences to the global reference sequence.
 *
 * @signature void join(stringSet, pos, joinConfig);
 *
 * @param[in,out] stringSet  The JournaledSet to join to.
 * @param[in]     pos        The position of the JournaledString within the string set.
 * @param[in]     joinConfig A @link JoinConfig @endlink object that specifies the method and the method's strategy
 *                           to compute the differences.
 */

// ----------------------------------------------------------------------------
// Function join()                                                [GlobalAlign]
// ----------------------------------------------------------------------------

template <typename TString, typename TPosition, typename TSpec>
inline void
join(StringSet<TString, Owner<JournaledSet> > & journalSet,
     TPosition journalIdx,
     JoinConfig<GlobalAlign<TSpec> > const & joinConfig)
{
   SEQAN_ASSERT_LT(journalIdx, static_cast<TPosition>(length(journalSet)));
   if (empty(host(journalSet)))
   {
       std::cerr << "No reference set! Join aborted!" << std::endl;
   }
   _joinInternal(value(journalSet, journalIdx), journalSet, joinConfig);
}

// ----------------------------------------------------------------------------
// Function join                                                  [GlobalChain]
// ----------------------------------------------------------------------------

template <typename TString, typename TPosition, typename TSpec>
inline void
join(StringSet<TString, Owner<JournaledSet> > & journalSet,
     TPosition journalIdx,
     JoinConfig<GlobalChain<TSpec> > const & joinConfig)
{
   SEQAN_ASSERT_LT(journalIdx, static_cast<TPosition>(length(journalSet)));

   if (empty(host(journalSet)))
   {
       std::cerr << "No reference set! Join aborted!" << std::endl;
   }
   _joinInternal(value(journalSet, journalIdx), journalSet, joinConfig);
}


// ----------------------------------------------------------------------------
// Function join                                                  [Simple Join]
// ----------------------------------------------------------------------------

template <typename TString, typename TPosition, typename TSpec>
inline void
join(StringSet<TString, Owner<JournaledSet> > & journalSet,
     TPosition journalIdx,
     TSpec const &)
{
   join(journalSet, journalIdx, JoinConfig<GlobalAlign<JournaledManhatten> >());
}

template <typename TString, typename TPosition>
inline void
join(StringSet<TString, Owner<JournaledSet> > & journalSet,
     TPosition journalIdx)
{
   join(journalSet, journalIdx, JoinConfig<GlobalAlign<JournaledManhatten> >());
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_JOIN_H_
