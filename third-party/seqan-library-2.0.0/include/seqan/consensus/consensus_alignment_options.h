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

#ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNMENT_OPTIONS_H_
#define INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNMENT_OPTIONS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConsensusAlignmentOptions
// ----------------------------------------------------------------------------

/*!
 * @class ConsensusAlignmentOptions
 * @headerfile <seqan/consensus.h>
 * @brief Configuration for @link consensusAlignment @endlink
 *
 * @signature struct ConsensusAlignmentOptions;
 */

struct ConsensusAlignmentOptions
{
    /*!
     * @var unsigned ConsensusAlignmentOptions::INVALID;
     * @brief Static member variable with a marker for invalid contig ids.
     */
    static const unsigned INVALID = (unsigned)-1;

    ConsensusAlignmentOptions() :
            contigID(INVALID), useContigID(true), usePositions(true), useGlobalAlignment(false),
            posDelta(30), runRealignment(true), verbosity(0), overlapMaxErrorRate(5), overlapMinLength(20),
            overlapMinCount(3), kMerSize(20), kMerMaxOcc(200),
            scoreMatch(5), scoreMismatch(-5), scoreGapOpen(-5), scoreGapExtend(-1)
    {}

    /*!
     * @var unsigned ConsensusAlignmentOptions::contigID;
     * @brief The id of the contig to compute the consensus for, defaults to <tt>INVALID</tt>.
     *
     * Set to <tt>INVALID</tt> for all.
     */
    unsigned contigID;

    /*!
     * @var bool ConsensusAlignmentOptions::useContigID;
     * @brief Whether or not to use the value of <tt>contigID</tt> (default is <tt>true</tt>).
     *
     * When this variable is set to <tt>false</tt> then the value of <tt>contigID</tt> is ignored and treated
     * as if it was <tt>INVALID</tt>.
     */
    bool useContigID;

    /*!
     * @var bool ConsensusAlignmentOptions::usePositions;
     * @brief Whether or not to use positions of the in the @link FragmentStore::alignedReadStore @endlink (default is <tt>true</tt>).
     *
     * When set to <tt>false</tt>, then the @link FragmentStore::alignedReadStore @endlink will be cleared and filled
     * with new entries, one for each read.
     */
    bool usePositions;

    bool useGlobalAlignment;

    /*!
     * @var unsigned ConsensusAlignmentOptions::posDelta;
     * @brief Positions are considered with an environment of <tt>posDelta</tt> (default is <tt>30</tt>).
     *
     * When positions are not used then this value is ignored.
     *
     * See @link consensusAlignment @endlink for more details.
    */
    unsigned posDelta;  // TODO(holtgrew): Rename to overlapWindowSize

    /*!
     * @var bool ConsensusAlignmentOptions::runRealignment;
     * @brief Perform a realignment using standard parameters, depending on the sequence length (default is <tt>true</tt>).
     *
     * Defaults to <tt>true</tt>.
     */
    bool runRealignment;

    int verbosity;

    int overlapMaxErrorRate;
    int overlapMinLength;
    int overlapMinCount;
    int kMerSize;
    int kMerMaxOcc;

    int scoreMatch;
    int scoreMismatch;
    int scoreGapOpen;
    int scoreGapExtend;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNMENT_OPTIONS_H_
