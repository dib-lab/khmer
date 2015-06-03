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

#ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_INTERFACE_H_
#define INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_INTERFACE_H_

#include <seqan/realign.h>
#include <seqan/store.h>

#include "consensus_alignment_options.h"

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
// Function consensusAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn consensusAlignment
 * @headerfile <seqan/consensus.h>
 * @brief Compute consensus alignment.
 *
 * @signature void consensusAlignment(store, options);
 *
 * @param[in,out] store   @link FragmentStore @endlink to use for consensus alignment computation.
 * @param[in]     options @link ConsensusAlignmentOptions @endlink with configuration.
 *
 * @throws ConsensusAlignerIllegalArgumentException in case of invalid arguments (e.g. two alignments for the same
 *                                                  read).
 *
 * This function computes a consensus alignment for a set of nucleic sequences that are stored in a @link FragmentStore
 * @endlink.  Often, consensus sequences are reads, but they could also be other sequences, such as RNA transcripts.
 * However, in the following description we call them reads.
 *
 * This function uses the @link FragmentStore::contigStore @endlink, @link FragmentStore::alignedReadStore @endlink,
 * and @link FragmentStore::readSeqStore @endlink members of <tt>store</tt>.
 *
 * Each read must have at most one entry in <tt>store.alignedReadStore</tt>.
 *
 * @section General Algorithm
 *
 * In the most common case, both contig ID information and position information is available.  The algorithm considers
 * all aligned reads on each contig.  For each read, all overlapping reads (with begin/end position extended
 * <tt>options.posDelta</tt> to the left/to the right) are considered and overlap alignments are computed.  This
 * pairwise alignment information is then used to compute a multiple sequence alignment (MSA).
 *
 * The resulting MSA is then refined by @link reAlignment @endlink (see @link ConsensusAlignmentOptions::runRealignment
 * @endlink).
 *
 * @section Using position information
 *
 * When position information is to be used then this will be used to generate fewer overlap alignments by only
 * considering possible overlaps in windows around each read alignment as described above.  Note that there can only be
 * at most one alignment for each read in the <tt>store.alignedReadStore</tt> and the end position must be greater than
 * or equal to the begin position, i.e., the alignment must be on the forward strand.
 *
 * Using position information also requires contig ID information.
 *
 * @section Using contigID information
 *
 * When contig ID information is to be used and position information is not to be used then the
 * <tt>consensusAlignment()</tt> will compute pairwise alignments between all pairs of reads on the same contig.
 *
 * When contigID information is to be used then only the reads with an entry in <tt>store.alignedReadStore</tt> are
 * considered.
 *
 * When contigID information is not used then an all-to-all pairwise alignment of all reads will be performed.
 *
 * <strong>Note</strong> that if reads having the same contig ID cannot be properly aligned and the MSA falls apart then
 * the reads for one "connected alignment component" are kept on the original contig while the rest are added to new
 * contigs that are appended to <tt>store.contigStore</tt>.
 *
 * @section Example
 *
 * The following example program takes a reference sequence and creates overlapping reads from it.  These are then added
 * to <tt>store</tt> with approximate positions (adding and subtracting a few positions).  The function
 * <tt>consensusAlignment()</tt> is then used to compute a MSA and the consensus sequence is stored in
 * <tt>store.contigStore[0].seq</tt>.
 *
 * @include demos/consensus/consensus_alignment.cpp
 *
 * The output is as follows:
 *
 * @include demos/consensus/consensus_alignment.cpp.stdout
 */

template <typename TSpec, typename TConfig>
void consensusAlignment(FragmentStore<TSpec, TConfig> & store,
                        ConsensusAlignmentOptions const & options)
{
    typedef FragmentStore<TSpec, TConfig> TFragmentStore;

    ConsensusAligner_<TFragmentStore> aligner(store, options);
    aligner.run();

    if (options.runRealignment)
        for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
            reAlignment(store, contigID, /*method=*/1, /*bandwidth=*/10, /*includeReference=*/false);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_ALIGNER_INTERFACE_H_
