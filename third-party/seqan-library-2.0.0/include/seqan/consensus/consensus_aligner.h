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

#ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_
#define INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_

#include <map>
#include <vector>
#include <stdexcept>

#include <seqan/store.h>

#include "consensus_alignment_options.h"
#include "consensus_builder.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConsensusAlignerIllegalArgumentException
// ----------------------------------------------------------------------------

/*!
 * @class ConsensusAlignerIllegalArgumentException
 * @headerfile <seqan/consensus.h>
 * @brief Thrown in ConsensusAlignerInputException on invalid arguments.
 *
 * @signature class ConsensusAlignerException;
 */

class ConsensusAlignerException : public std::runtime_error
{
public:
    explicit ConsensusAlignerException(std::string const & whatArg) : std::runtime_error(whatArg)
    {}
};

// ----------------------------------------------------------------------------
// Class ConsensusAligner_
// ----------------------------------------------------------------------------

// Implementation helper of the consensusAlignment() function.

template <typename TFragmentStore>
class ConsensusAligner_
{
public:
    ConsensusAligner_(TFragmentStore & store, ConsensusAlignmentOptions const & options) :
            store(store), options(options)
    {}

    void run();

private:

    // Checks that there is at most one alignment per read in the store.
    void checkAlignmentMultiplicity() const;
    // Normalize reverse-complemented alignments and store the read IDs in rcIDs.
    void normalizeRCReads();
    // Restore read orientation of RC reads again.
    void restoreRCReads();

    // The ids of the reads that had a reverse-complemented alignment.
    std::set<unsigned> rcIDs;
    // The FragmentStore to use for consensus computation.
    TFragmentStore & store;
    // The configuration of the consensus alignment.
    ConsensusAlignmentOptions const & options;
};

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::normalizeRCReads()
{
    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
        if (store.alignedReadStore[i].beginPos > store.alignedReadStore[i].endPos)
        {
            unsigned readID = store.alignedReadStore[i].readId;
            rcIDs.insert(readID);
            reverseComplement(store.readSeqStore[readID]);
            std::swap(store.alignedReadStore[i].beginPos, store.alignedReadStore[i].endPos);
        }
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::restoreRCReads()
{
    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
    {
        unsigned readID = store.alignedReadStore[i].readId;
        if (rcIDs.count(readID))
        {
            std::swap(store.alignedReadStore[i].beginPos, store.alignedReadStore[i].endPos);
            reverseComplement(store.readSeqStore[readID]);
        }
    }
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::checkAlignmentMultiplicity() const
{
    std::set<unsigned> seen;
    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
        if (seen.count(store.alignedReadStore[i].readId))
            throw ConsensusAlignerException("Read with id has more than one alignment.");
        else
            seen.insert(store.alignedReadStore[i].readId);
}

template <typename TFragmentStore>
inline void ConsensusAligner_<TFragmentStore>::run()
{
    // Check that there only is one alignment per read and escape via exception otherwise.
    checkAlignmentMultiplicity();
    // Normalize reads with RC alignment.
    normalizeRCReads();

    if (options.verbosity >= 1)
        std::cerr << "computing overlap infos...\n";
    // Build overlap infos, based on position and contig ID if configured to do so.
    std::vector<OverlapInfo_> overlapInfos;
    OverlapInfoComputation_<TFragmentStore> ovlInfoHelper(store, options);
    ovlInfoHelper.run(overlapInfos);
    if (options.verbosity >= 1)
        std::cerr << "=> " << overlapInfos.size() << " overlap infos\n";

    // TODO(holtgrew): Pick pairwise best alignments.

    if (options.verbosity >= 1)
        std::cerr << "building consensus...\n";
    // Consensus builder.
    ConsensusBuilder_<TFragmentStore> consBuilder(options);
    consBuilder.run(store, overlapInfos);

    // Restore RC'ed reads.
    restoreRCReads();
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_H_
