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

#ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_BUILDER_H_
#define INCLUDE_SEQAN_CONSENSUS_CONSENSUS_BUILDER_H_

#include <map>

#include <seqan/graph_types.h>
#include <seqan/store.h>

#include "consensus_alignment_options.h"
#include "overlapper.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class ConsensusBuilder_
// --------------------------------------------------------------------------

template <typename TFragmentStore>
class ConsensusBuilder_
{
public:
    ConsensusBuilder_(ConsensusAlignmentOptions const & options) : options(options)
    {}

    void run(TFragmentStore & store, std::vector<OverlapInfo_> const & overlapInfos) const;

private:

    // Glorified pair<unsigned, int> with more verbose name for readID and offset.
    struct AlignmentTarget
    {
        AlignmentTarget() : readID((unsigned)-1), offset(0) {}
        AlignmentTarget(unsigned readID, int offset) : readID(readID), offset(offset) {}

        unsigned readID;
        int offset;
    };

    // Glorified triple for overlap computation.
    struct PositionedRead
    {
        PositionedRead() : beginPos(0), endPos(0), readID((unsigned)-1) {}
        PositionedRead(int beginPos, int endPos, unsigned readID) :
            beginPos(beginPos), endPos(endPos), readID(readID)
        {}

        int beginPos;
        int endPos;
        unsigned readID;

        bool operator<(PositionedRead const & other) const
        {
            return (std::make_pair(std::make_pair(beginPos, endPos), readID) <
                    std::make_pair(std::make_pair(other.beginPos, other.endPos), other.readID));
        }
    };

    // Shortcut for fragments.
    typedef String<Fragment<> > TFragments;

    // Compute pairwise overlaps for the reads of each contig in cg.
    //
    // In case that there are regions higher than options.startConsensusCompression, we add subsample these reads and
    // write out an entry (other read id, offset) to use as a position marker in the resulting MSA.
    void computeOverlaps(
            TFragments & fragments,
            String<int> & scores,
            Graph<Undirected<double> > & distances,
            TFragmentStore const & store,
            std::vector<OverlapInfo_> const & overlapInfos) const;

    // Configuration to use for consensus building.
    ConsensusAlignmentOptions options;
};

template <typename TFragmentStore>
void ConsensusBuilder_<TFragmentStore>::computeOverlaps(
            typename ConsensusBuilder_<TFragmentStore>::TFragments & fragments,
            String<int> & scores,
            Graph<Undirected<double> > & distances,
            TFragmentStore const & store,
            std::vector<OverlapInfo_> const & overlapInfos) const
{
    // Configure overlapper, no limit in error rate or minimum length.
    OverlapperOptions_ ovlOptions;
    ovlOptions.logging = (options.verbosity >= 3);
    ovlOptions.overlapErrorRate = options.overlapMaxErrorRate;
    ovlOptions.overlapMinLength = options.overlapMinLength;  // accept all overlaps
    typedef typename Value<typename TFragmentStore::TReadSeqStore const>::Type TReadSeq;
    Overlapper_<String<Fragment<> >, TReadSeq> overlapper(ovlOptions);

    for (std::vector<OverlapInfo_>::const_iterator it = overlapInfos.begin(); it != overlapInfos.end(); ++it)
    {
        Overlap_ ovl;
        TFragments frags;
        if (overlapper.computeOverlap(ovl, frags, store.readSeqStore[it->seq0], store.readSeqStore[it->seq1],
                                      OverlapCandidate_(it->seq0, it->seq1, it->pos1 - it->numErrors,
                                                        it->pos1 + it->numErrors)))
        {
            append(fragments, frags);
            int ovlScore = ovl.length() - ovl.errors;  // pseudo-score
            resize(scores, length(fragments), ovlScore);
            int quality = (ovl.length() > 0) ? (ovl.length() - ovl.errors) * 100 / ovl.length() : 0;
            addEdge(distances, ovl.seq0, ovl.seq1, quality);
        }
    }
}

// ---------------------------------------------------------------------------
// Function alignmentGraphToSmoothFragmentStore()
// ---------------------------------------------------------------------------

// Note that this function relies on the "all mated, adjacent reads" assumption.

template <typename TFragmentStore, typename TSequence, typename TCargo, typename TSetSpec, typename TSpec>
bool alignmentGraphToFragmentStore(TFragmentStore & store,
                                   Graph<Alignment<StringSet<TSequence, TSetSpec>, TCargo, TSpec> > const & g,
                                   Graph<Undirected<double> > const & distances,
                                   String<unsigned> const & component,
                                   String<unsigned> const & order,
                                   unsigned numComponents,
                                   bool logging)
{
    bool const DEBUG_INCONSISTENT_LEN = false;

    // std::cerr << ">>>>>>>>>>>>\n<<<<<<<<<<<<<<<<\n";
    // NOTE: seqToCluster is indexed by POSITION in the read set of g and not by the ID.

    // TODO(holtgrew): This function is very similar to the updateStoreFromAlignmentGraph computeProfiles functions. Maybe we can share the commonality?
    typedef Graph<Alignment<StringSet<TSequence, TSetSpec>, TCargo, TSpec> > TAlignmentGraph;

    // -----------------------------------------------------------------------
    // Get connected components of distances / read alignment clusters.
    // -----------------------------------------------------------------------

    // Each cluster corresponds to a contig.

    // A cluster is a CC in the graph where each sequences is a vertex and two vertices are connected if they have an
    // overlap alignment.
    String<unsigned> seqToCluster;
    if (logging)
        std::cerr << "# vertices: " << numVertices(distances) << "\n"
                  << "# edges: " << numEdges(distances) << "\n";
    unsigned numClusters = connectedComponents(seqToCluster, distances);
    if (logging)
        std::cerr << "# clusters: " << numClusters << std::endl
                  << "# components: " << numComponents << std::endl;
    clear(store.contigStore);
    resize(store.contigStore, numClusters);
    String<unsigned> contigLengths;
    resize(contigLengths, numClusters, 0);

    for (unsigned i = 0; i < numClusters; ++i)
    {
        std::stringstream ss;
        ss << "contig_" << i;
        appendValue(store.contigNameStore, ss.str());
    }

    // -----------------------------------------------------------------------
    // Visit components in topological order and generate profile sequences.
    // -----------------------------------------------------------------------

    // Get mapping from component to vertices.
    String<String<unsigned> > componentVertices;
    resize(componentVertices, numComponents);
    typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIterator;
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV))
        appendValue(componentVertices[getProperty(component, *itV)], *itV);

    // For each cluster, the currently overlapping reads.
    std::vector<std::set<unsigned> > activeReads(numClusters);

    std::vector<unsigned> gapCount(length(stringSet(g)), 0);

    // Recreate alignedReadStore.
    clear(store.alignedReadStore);
    resize(store.alignedReadStore, length(store.readSeqStore));

    // Iterate vertices in topological order.
    for (typename Iterator<String<unsigned> const, Rooted>::Type it = begin(order, Rooted()); !atEnd(it); goNext(it))
    {
        unsigned c = *it;     // Current component.
        unsigned fLen = fragmentLength(g, front(componentVertices[c]));
        for (unsigned i = 1; i < length(componentVertices[c]); ++i)
            SEQAN_ASSERT_EQ(fragmentLength(g, front(componentVertices[c][0])),
                            fragmentLength(g, front(componentVertices[c][i])));
        unsigned cl = seqToCluster[idToPosition(stringSet(g), sequenceId(g, front(componentVertices[c])))];  // Current cluster/contig.

        // Update contig lengths.
        unsigned from = contigLengths[cl];
        contigLengths[cl] += fLen;
        if (DEBUG_INCONSISTENT_LEN)
            std::cerr << "==== c == " << c << "\n" << "     from == " << from << "\n";

        // The currently active reads that we see in this round.  Required for inserting gaps below.
        std::set<unsigned> seen;
        std::set<unsigned> done;

        // Insert gaps.
        typedef typename Iterator<String<unsigned>, Rooted>::Type TDescIt;
        for (TDescIt itV = begin(componentVertices[c], Rooted()); !atEnd(itV); goNext(itV))
        {
            unsigned idx = idToPosition(stringSet(g), sequenceId(g, *itV));
            seen.insert(idx);
            unsigned fBeg = fragmentBegin(g, *itV);
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "           fBeg of " << *itV << " is " << fBeg << " (idx == " << idx << ")\n";

            // Register sequence as supporting in profile cl starting at position from in profile.
            if (fBeg == 0u)
            {
                activeReads[cl].insert(idx);
                store.alignedReadStore[idx].id = idx;
                store.alignedReadStore[idx].readId = idx;
                store.alignedReadStore[idx].contigId = cl;
                store.alignedReadStore[idx].beginPos = from;
                store.alignedReadStore[idx].endPos = from;
                // store.alignedReadStore[idx].pairMatchId = idx / 2;  // TODO(holtgrew): Set from read pair info.
                if (DEBUG_INCONSISTENT_LEN)
                {
                    std::cerr << "store.alignedReadStore[" << idx << "].beginPos == " << from << " | *itV == " << *itV << "\n";
                    std::cerr << "store.alignedReadStore[" << idx << "].endPos (= endPos) == " << from << "| *itV == " << *itV << "\n";
                }
            }
            store.alignedReadStore[idx].endPos = from + fLen;
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "store.alignedReadStore[" << idx << "].endPos = " << from << " + " << fLen << " == " << (from + fLen) << "| *itV == " << *itV << "\n";

            unsigned fEnd = fBeg + fLen;

            if (fEnd == length(stringSet(g)[idx]))
                done.insert(idx);
        }

        // Get not seen reads.
        typedef std::set<unsigned>::iterator TSetIt;
        std::set<unsigned> notSeen;
        for (TSetIt it = activeReads[cl].begin(); it != activeReads[cl].end(); ++it)
            notSeen.insert(*it);
        for (TSetIt it = seen.begin(); it != seen.end(); ++it)
            notSeen.erase(*it);
        // Insert gaps into these reads.
        for (TSetIt itS = notSeen.begin(); itS != notSeen.end(); ++itS)
        {
            typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
            typedef typename Value<TAlignedReadStore>::Type TAlignedRead;
            typedef typename TAlignedRead::TGapAnchors TGapAnchors;
            typedef typename TFragmentStore::TReadSeq TReadSeq;
            SEQAN_ASSERT_NOT(empty(store.readSeqStore[*itS]));
            Gaps<TReadSeq, AnchorGaps<TGapAnchors> > gaps(store.readSeqStore[*itS], store.alignedReadStore[*itS].gaps);
            insertGaps(gaps, from - store.alignedReadStore[*itS].beginPos, fLen);
            store.alignedReadStore[*itS].endPos += fLen;
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "store.alignedReadStore[" << *itS << "].endPos += " << fLen << " == " << store.alignedReadStore[*itS].endPos << "\n";
            gapCount[*itS] += fLen;
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "gapCount[" << *itS << "] == " << gapCount[*itS] << "\n";
        }

        // Deactive done reads.
        for (TSetIt it = done.begin(); it != done.end(); ++it)
            activeReads[cl].erase(*it);
    }

// #if SEQAN_ENABLE_DEBUG
    {
        // Check for consistency.
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
        typedef typename TFragmentStore::TReadSeq           TReadSeq;

        TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
        for (TAlignedReadIter it2 = begin(store.alignedReadStore, Standard()); it2 != itEnd; ++it2)
        {
            typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
            TReadGaps readGaps(store.readSeqStore[it2->readId], it2->gaps);
            SEQAN_ASSERT_EQ(length(readGaps) - length(store.readSeqStore[it2->readId]), gapCount[it2->readId]);
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "READ GAPS\t" << (it2 - begin(store.alignedReadStore, Standard())) << "\t>>>" << readGaps << "<<< (" << length(readGaps) << ")\n"
                          << "  beginPos == " << it2->beginPos << ", endPos == " << it2->endPos << ", gapCount == " << gapCount[it2->readId] << "\n";
            if ((unsigned)std::abs(it2->endPos - it2->beginPos) != length(readGaps))
            {
                SEQAN_FAIL("Inconsistent begin/endPos");
            }
        }
    }
// #endif  // #if SEQAN_ENABLE_DEBUG

    return true;
}

template <typename TFragmentStore, typename TSequence, typename TCargo, typename TSetSpec, typename TSpec>
bool alignmentGraphToFragmentStore(TFragmentStore & store,
                                   Graph<Alignment<StringSet<TSequence, TSetSpec>, TCargo, TSpec> > const & g,
                                   Graph<Undirected<double> > const & distances,
                                   bool logging)
{
    typedef std::map<unsigned, unsigned> TComponentLength;

    // -----------------------------------------------------------------------
    // Compute connected components and get topological sorting of them.
    // -----------------------------------------------------------------------
      String<unsigned> component;
      String<unsigned> order;
      TComponentLength componentLength;
    if (empty(g))
        return true;  // Nothing to do for empty graphs.
      if (!convertAlignment(g, component, order, componentLength))
        return false;
    unsigned numComponents = length(order);

    return alignmentGraphToFragmentStore(store, g, distances, component, order, numComponents, logging);
}

template <typename TFragmentStore>
void ConsensusBuilder_<TFragmentStore>::run(TFragmentStore & store,
                                            std::vector<OverlapInfo_> const & overlapInfos) const
{
    // The alignments as required by consensus MSA module.
    TFragments fragments;
    String<int> scores;
    Graph<Undirected<double> > distances;
    _resizeWithRespectToDistance(distances, length(store.readSeqStore));

    // Obtain alignments according to contig graph.
    computeOverlaps(fragments, scores, distances, store, overlapInfos);

    // Build copy of dependent reads.
    typedef typename Value<typename TFragmentStore::TReadSeqStore>::Type TReadSeq;
    StringSet<TReadSeq> seqs;
    for (unsigned i = 0; i < length(store.readSeqStore); ++i)
        appendValue(seqs, store.readSeqStore[i]);

    // Build the alignment graph.
    Score<int, Simple> msaScoringScheme(2, -6, -4, -9);  // TODO(holtgrew): Get from options!
    typedef StringSet<TReadSeq, Dependent<> > TDepReadSet;
    typedef Graph<Alignment<TDepReadSet, unsigned> > TInGraph;
    TDepReadSet depSeqs(seqs);
    TInGraph inGraph(depSeqs);
    buildAlignmentGraph(fragments, scores, inGraph, msaScoringScheme, ReScore());

    if (options.verbosity >= 2)
        std::cerr << "Building alignment graph\n"
                  << "  # fragments: " << length(fragments) << "\n"
                  << "  # seqs: " << length(seqs) << "\n";

    // Perform triplet library extension.
    if (length(seqs) < 1000/*(unsigned)options.texStopCount*/)  // TODO(holtgrew): Clean up.
        tripletLibraryExtension(inGraph);

    Graph<Tree<double> > guideTree;
    Graph<Undirected<double> > dCopy(distances);
    upgmaTree(dCopy, guideTree);

    // Perform progressive alignment.
    if (options.verbosity >= 2)
        std::cerr << "progressive alignment...\n";
    TInGraph graph(depSeqs);
    assignStringSet(graph, stringSet(inGraph));
    progressiveAlignment(inGraph, guideTree, graph);

    // Build alignment graph.
    if (options.verbosity >= 2)
        std::cerr << "build AG..\n";

    String<unsigned> component;
    String<unsigned> order;
    std::map<unsigned, unsigned> componentLength;
    SEQAN_ASSERT_NOT(empty(graph));
    bool b = convertAlignment(graph, component, order, componentLength);
    (void) b;
    SEQAN_ASSERT(b);
    unsigned numComponents = length(order);
    if (options.verbosity >= 2)
        std::cerr << "AG\n"
                  << "  numVertices = " << numVertices(graph) << "\n"
                  << "  numEdges = " << numEdges(graph) << "\n"
                  << "  componentLength.size() = " << componentLength.size() << "\n"
                  << "  numComponents = " << numComponents << "\n";

    if (options.verbosity >= 2)
        std::cerr << "build store..\n";
    alignmentGraphToFragmentStore(store, graph, distances, component, order, numComponents,
                                  /*logging=*/(options.verbosity >= 3));

    if (options.verbosity >= 2)
    {
        std::cerr << "CONSENSUS RESULT\n";
        typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignedReadIter;
        for (TAlignedReadIter it = begin(store.alignedReadStore, Standard()); it != end(store.alignedReadStore, Standard()); ++it)
            std::cerr << "contigID=" << it->contigId
                      << ", beginPos=" << std::min(it->beginPos, it->endPos)
                      << ", readID=" << it->readId
                      << ", seq=" << store.readSeqStore[it->readId]
                      << "\n";
        std::cerr << "\n";


    }

    if (options.verbosity >= 2)
    {
        std::cerr << "before realignment\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
        {
            int endPos = 0;
            for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
                if (store.alignedReadStore[i].contigId == contigID)
                    endPos = std::max(endPos, (int)store.alignedReadStore[i].endPos);
            std::cerr << ">contig_" << contigID << "\n";
            printAlignment(std::cerr, layout, store, /*contigID=*/contigID, /*beginPos=*/0, /*endPos=*/endPos, 0, 30);
        }
    }

    if (options.verbosity >= 2)
        std::cerr << "realigning...\n";
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        reAlignment(store, i, 1, 10, false);


    if (options.verbosity >= 2)
    {
        std::cerr << "after realignment\n";
        seqan::AlignedReadLayout layout;
        layoutAlignment(layout, store);
        for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
        {
            int endPos = 0;
            for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
                if (store.alignedReadStore[i].contigId == contigID)
                    endPos = std::max(endPos, (int)store.alignedReadStore[i].endPos);
            std::cerr << ">contig_" << contigID << "\n";
            printAlignment(std::cerr, layout, store, /*contigID=*/contigID, /*beginPos=*/0, /*endPos=*/endPos, 0, 30);
        }
    }
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_CONSENSUS_CONSENSUS_BUILDER_H_
