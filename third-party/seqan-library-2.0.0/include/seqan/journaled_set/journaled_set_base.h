// ==========================================================================
//                               journaled_set
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

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_BASE_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Spec JournaledSet
// ----------------------------------------------------------------------------

/*!
 * @class JournaledSet
 * @extends StringSet
 * @headerfile <seqan/journaled_set.h>
 *
 * @signature template <typename TString>
 *            class StringSet<TString, Owner<JournaledSet> >;
 * @tparam TString The string type. Types: @link String @endlink, @link JournaledString @endlink
 *
 * @brief A @link StringSet @endlink storing the strings as members.  It can store a global reference sequence to which all members can
 * be journaled if they are of type @link JournaledString @endlink.
 *
 * The strings are internally stored in a <tt>String&lt;TString&gt;</tt> object and the character position type is a @link
 * Pair @endlink <tt>(seqNo, seqOfs)</tt> where seqNo identifies the string within the string set and seqOfs identifies
 * the position within this string.
 *
 * The global reference is of type <tt>Host&lt;TString&gt;</tt>. Only strings of type @link JournaledString @endlink or
 * <tt>Host&lt;</tt>@link JournaledString @endlink<tt>&gt;</tt> can be used for the advanced functionality supported by
 * this string set.
 */

struct JournaledSet_;
typedef Tag<JournaledSet_> JournaledSet;


// ----------------------------------------------------------------------------
// Class JournalTraceBuffer
// ----------------------------------------------------------------------------

template <typename TString>
class JournalTraceBuffer;

/*!
 * @defgroup JoinStrategiesTags Join Strategies Tags
 * @brief Tags used for selecting journaling strategies when joining a JournaledString to a global reference sequence.
 *
 *
 * @tag JoinStrategiesTags#JournaledManhattan
 * @headerfile <seqan/sequence_journaled.h>
 * @brief Constructs a @link JournaledString @endlink based on Manhattan distance.
 *
 * This strategy is very fast on the cost of memory.
 *
 * @signature typedef Tag<JournaledManhattan_> JournaledManhattan.
 *
 *
 * @tag JoinStrategiesTags#JournaledCompact
 * @headerfile <seqan/sequence_journaled.h>
 * @brief Computes an optimal alignment to construct a @link JournaledString @endlink.
 *
 * This strategy is slow but depending on the scoring function minimizes the memory requirements for the computed @link
 * JournaledSet @endlink.
 *
 * @signature typedef Tag<JournaledCompact_> JournaledCompact;
 */

// ----------------------------------------------------------------------------
// Tag JournaledManhatten
// ----------------------------------------------------------------------------

struct JournaledManhatten_;
typedef Tag<JournaledManhatten_> JournaledManhatten;

// ----------------------------------------------------------------------------
// Tag JournaledCompact
// ----------------------------------------------------------------------------

struct JournaledCompact_;
typedef Tag<JournaledCompact_> JournaledCompact;


// ----------------------------------------------------------------------------
// Spec GlobalAlign
// ----------------------------------------------------------------------------

/*!
 * @class GlobalAlign
 * @extends JoinConfig
 *
 * @headerfile <seqan/journaled_set.h>
 * @brief Selects a global alignment method to join a @link JournaledString @endlink to a global reference sequence.
 *
 * @signature template <[typename TStrategy]>
 *            struct GlobalAlign;
 *
 * @tparam TStrategy The strategy used to compute the journal (can be one of @link JoinStrategiesTags#JournaledCompact
 *                   @endlink and @link JoinStrategiesTags#JournaledManhattan @endlink, defaults to @link
 *                   JoinStrategiesTags#JournaledManhattan @endlink).
 *
 * If @link JoinStrategiesTags#JournaledManhattan @endlink is selected, then the resulting @link JournaledString
 * @endlink consists of one insertion node covering the complete joined sequence.
 */

template <typename TSpec = JournaledManhatten>
struct GlobalAlign{};

// ----------------------------------------------------------------------------
// Spec GlobalChain
// ----------------------------------------------------------------------------

/*!
 * @class GlobalChain
 * @extends JoinConfig
 *
 * @headerfile <seqan/journaled_set.h>
 * @brief Selects an anchor-based method to join a @link JournaledString @endlink to a global reference sequence.
 *
 * @signature template <[typename TStrategy]>
 *            struct GlobalChain;
 *
 * @tparam TStrategy The strategy used to compute the journal (can be one of @link JoinStrategiesTags#JournaledCompact
 *                   @endlink and @link JoinStrategiesTags#JournaledManhattan @endlink, defaults to @link
 *                   JoinStrategiesTags#JournaledManhattan @endlink).
 *
 * The @link JoinStrategiesTags#JournaledManhattan @endlink strategy fills the gaps between the anchors with a single
 * insertion node whil the corresponding part of the reference sequence is deleted.
 */

template <typename TSpec = JournaledManhatten>
struct GlobalChain{};

// ----------------------------------------------------------------------------
// Class JoinConfig
// ----------------------------------------------------------------------------

/*!
 * @class JoinConfig
 * @headerfile <seqan/journaled_set.h>
 * @brief Specifies the strategy and all necessary parameters used to journal a sequence to a reference sequence.
 *
 * @signature template <[typename TMethod]>
 *            struct JoinConfig;
 *
 * @tparam TMethod The method type. Types: @link GlobalAlign @endlink, @link GlobalChain @endlink
 *
 * SeqAn offers two general methods to compute the journal.  The first one uses a @link globalAlignment @endlink
 * function and the second one uses an anchor based approach.
 */

template <typename TSpec = GlobalAlign<JournaledManhatten> >
struct JoinConfig{};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


}  // namespace seqan

#endif  // INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_BASE_H_
