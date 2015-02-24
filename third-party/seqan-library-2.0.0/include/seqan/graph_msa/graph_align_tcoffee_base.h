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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Distance unity
//////////////////////////////////////////////////////////////////////////////

//static const int SEQAN_DISTANCE_UNITY = 1 << 20;
static const int SEQAN_DISTANCE_UNITY = 1;


//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*!
 * @defgroup AlignmentGraphCombinationTags Alignment Graph Combination
 * @brief A tag to specify how to combine @link AlignmentGraph AlignmentGraphs @endlink.
 *
 *
 * @tag AlignmentGraphCombinationTags#FractionalScore
 * @headerfile <seqan/graph_msa.h>
 * @brief Rescore matches with appropriate fractional scores.
 *
 * @signature typedef Tag<FractionalScore_> const FractionalScore;
 *
 *
 * @tag AlignmentGraphCombinationTags#FrequencingCounting
 * @headerfile <seqan/graph_msa.h>
 * @brief Rescores matches with the frequency count for this edge.
 *
 * @signature typedef Tag<FrequencingCounting_> const FrequencingCounting;
 *
 *
 * @tag AlignmentGraphCombinationTags#ReScore
 * @headerfile <seqan/graph_msa.h>
 * @brief Rescore the matches after segment-match refinement.
 *
 * @signature typedef Tag<FractionalScore_> const ReScore;
 */

struct FractionalScore_;
typedef Tag<FractionalScore_> const FractionalScore;


struct FrequencyCounting_;
typedef Tag<FrequencyCounting_> const FrequencyCounting;

struct ReScore_;
typedef Tag<ReScore_> const ReScore;


//////////////////////////////////////////////////////////////////////////////
// Generating an alignment graph from segment matches
//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn buildAlignmentGraph
 * @headerfile <seqan/graph_msa.h>
 * @brief Builds an @link AlignmentGraph @endlink from a seqt of input alignments.

 * @signature void buildAlignmentGraph(matches[, score], outGraph, tag);
 * @signature void buildAlignmentGraph(matches, score, outGraph, scoreType, ReScore);
 *
 * @param[in]  matches   A @link String @endlink of alignments, e.g. @link Fragment @endlink or @link Align @endlink.
 * @param[in]  score     A @link String @endlink of score values.
 * @param[out] outGraph  An @link AlignmentGraph @endlink with a previously filled @link DependentStringSet @endlink.
 * @param[in]  scoreType A @link Score @endlink type.
 * @param[in]  tag       A tag indicating the score of matches.  See @link AlignmentGraphCombinationTags @endlink.
 *
 * Calls the function @link matchRefinement @endlink and adapts the scores according to <tt>tag</tt>.
 */


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TCargo, typename TSpec>
inline void
buildAlignmentGraph(String<TFragment, TSpec1>& matches,
                    String<TScoreValue, TSpec2>& scores,
                    Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
                    FractionalScore)
{
    SEQAN_CHECKPOINT
    typedef String<TFragment, TSpec1> TFragmentString;
    typedef typename Iterator<TFragmentString, Standard>::Type TFragmentStringIter;
    typedef String<TScoreValue, TSpec2> TScoreValues;
    typedef typename Iterator<TScoreValues, Standard>::Type TScoreValuesIter;
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
    typedef typename Size<TFragmentString>::Type TSize;
    typedef typename Id<TOutGraph>::Type TId;
    //typedef typename EdgeDescriptor<TOutGraph>::Type TEdgeDescriptor;
    typedef typename VertexDescriptor<TOutGraph>::Type TVertexDescriptor;

    // Initialization
    clearVertices(outGraph);
    TStringSet& strSet = stringSet(outGraph);

    // Segment-match refinement
    matchRefinement(matches,strSet,outGraph);

    // Clear edge-weights
    typedef typename Iterator<TOutGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itE(outGraph);
    for(;!atEnd(itE);goNext(itE)) cargo(value(itE)) = 1;

    // Adapt the scores
    TFragmentStringIter it = begin(matches, Standard() );
    TFragmentStringIter endIt = end(matches, Standard() );
    TScoreValuesIter scoreIt = begin(scores, Standard() );
    for(; it != endIt; ++it, ++scoreIt) {
        TId id1 = sequenceId(*it,0);
        TId id2 = sequenceId(*it,1);
        TSize pos1 = fragmentBegin(*it, id1);
        TSize pos2 = fragmentBegin(*it, id2);
        TSize fragLen = fragmentLength(*it, id1);
        TSize end1 = pos1 + fragLen;
        while(pos1 < end1) {
            TVertexDescriptor p1 = findVertex(outGraph, id1, pos1);
            TVertexDescriptor p2 = findVertex(outGraph, id2, pos2);
            TSize vertexLen = fragmentLength(outGraph, p1);
            cargo(findEdge(outGraph, p1, p2)) += (TCargo) ((vertexLen * (*scoreIt)) / fragLen);
            pos1 += vertexLen;
            pos2 += vertexLen;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TStringSet, typename TCargo, typename TSpec>
inline void
buildAlignmentGraph(String<TFragment, TSpec1>& matches,
                    Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
                    FrequencyCounting)
{
    //typedef String<TFragment, TSpec1> TFragmentString;
    //typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
    //typedef typename Size<TFragmentString>::Type TSize;

    // Initialization
    clearVertices(outGraph);
    TStringSet& strSet = stringSet(outGraph);

    // Segment-match refinement
    matchRefinement(matches,strSet,outGraph);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TCargo, typename TSpec>
inline void
buildAlignmentGraph(String<TFragment, TSpec1>& matches,
                    String<TScoreValue, TSpec2>&,
                    Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
                    FrequencyCounting)
{
    buildAlignmentGraph(matches, outGraph, FrequencyCounting() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TScoreType, typename TSize, typename TSpec2, typename TScoreString, typename TScoreValue>
inline void
_scoreMatches(StringSet<TString, TSpec> const & seqSet,
              TScoreType const & scType,
              String<Fragment<TSize, ExactFragment<> >, TSpec2> const & matches,
              TScoreString & scores,
              TScoreValue offset)
{
    typedef String<Fragment<TSize, ExactFragment<> >, TSpec2> const TFragmentString;
    typedef typename Id<typename Value<TFragmentString>::Type>::Type TId;
    typedef typename Iterator<TFragmentString, Standard>::Type TFragmentStringIter;
    typedef typename Iterator<TString const, Standard>::Type TStringIter;
    typedef typename Iterator<TScoreString, Standard>::Type TScoreStringIter;
    resize(scores, length(matches));

    // Get the scores
    TFragmentStringIter itF = begin(matches, Standard() );
    TFragmentStringIter itFEnd = end(matches, Standard() );
    TScoreStringIter itSc = begin(scores, Standard() );
    TId id1 = 0; TId id2 = 0;
    TSize pos1 = 0; TSize pos2 = 0;    TSize fragLen = 0;
    for(; itF != itFEnd; ++itF, ++itSc) {
        id1 = sequenceId(*itF,0);
        id2 = sequenceId(*itF,1);
        pos1 = fragmentBegin(*itF, id1);
        pos2 = fragmentBegin(*itF, id2);
        fragLen = fragmentLength(*itF, id1);
        TStringIter itS1 = begin(seqSet[idToPosition(seqSet, id1)], Standard() );
        itS1 += pos1;
        TStringIter itS2 = begin(seqSet[idToPosition(seqSet, id2)], Standard() );
        itS2 += pos2;
        *itSc = 0;
        for(TSize i = 0; i<fragLen; ++i, ++itS1, ++itS2)
            *itSc += offset + score(scType, *itS1, *itS2);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TScoreType, typename TFragment, typename TSpec2, typename TScoreString>
inline void
_scoreMatches(StringSet<TString, TSpec> const& seqSet,
              TScoreType const& scType,
              String<TFragment, TSpec2> const& matches,
              TScoreString& scores)
{
    SEQAN_CHECKPOINT
    _scoreMatches(seqSet, scType, matches, scores, (typename Value<TScoreString>::Type) 10);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void
buildAlignmentGraph(String<TFragment, TSpec1>& matches,
                    String<TScoreValue, TSpec2>& scores,
                    Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
                    TScore const& scType,
                    ReScore)
{
    // ReScore
    _scoreMatches(stringSet(outGraph), scType, matches, scores);

    // Use fractinal score
    buildAlignmentGraph(matches, scores, outGraph, FractionalScore() );
}


//////////////////////////////////////////////////////////////////////////////
// Consistency: Triplet extension
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



template<typename TVertexDescriptor, typename TCargo>
struct MsaEdgeCargo_ {
 public:
     TVertexDescriptor v1;
     TVertexDescriptor v2;
     TCargo c;

     MsaEdgeCargo_() {}


     MsaEdgeCargo_(TVertexDescriptor vert1, TVertexDescriptor vert2, TCargo carg) :
     v1(vert1), v2(vert2), c(carg) {}
};

//////////////////////////////////////////////////////////////////////////////

template<typename TVertexDescriptor, typename TCargo>
struct LessMsaEdgeCargo_ :
    public std::binary_function<TVertexDescriptor, TCargo, bool>
{
    inline bool
    operator() (MsaEdgeCargo_<TVertexDescriptor, TCargo> const& a1,
                MsaEdgeCargo_<TVertexDescriptor, TCargo> const& a2) const
    {
        return (a1.v1 == a2.v1) ? (a1.v2 < a2.v2) : (a1.v1 < a2.v1);
    }
};

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn tripletLibraryExtension
 * @headerfile <seqan/graph_msa.h>
 * @brief Performs a full or group-based consistency extension.
 *
 * @signature void tripletLibraryExtension(graph[, guideTree, minMembers]);
 *
 * @param[in,out] graph      An @link AlignmentGraph @endlink.
 * @param[in]     guideTree  The guide tree to use.
 * @param[in]     minMembers Minimum number of sequences per group.  If a guide tree and a minimum number of members
 *                           if given then the triplet extension is limited to groups of sequences.
 *
 * The running time is quadratic in the number of pairwise edges.
 */

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
tripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    //typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

    // Store all edges
    typedef std::pair<TVertexDescriptor, TVertexDescriptor> TNewEdge;
    typedef std::map<TNewEdge, TCargo> TEdgeMap;
    typedef typename TEdgeMap::iterator TEdgeMapIter;
    TEdgeMap newEMap;
    typedef MsaEdgeCargo_<TVertexDescriptor, TCargo> TEdge;
    typedef String<TEdge> TEdgeString;
    TEdgeString fullEdges;
    TEdgeIterator itE(g);
    for(;!atEnd(itE);goNext(itE)) {
        TVertexDescriptor sV = sourceVertex(itE);
        TVertexDescriptor tV = targetVertex(itE);
        TCargo c = cargo(*itE);
        appendValue(fullEdges, TEdge(sV, tV, c), Generous());
        appendValue(fullEdges, TEdge(tV, sV, c), Generous());
        newEMap.insert(std::make_pair(TNewEdge(sV, tV), c));
    }
    clearEdges(g);
    std::sort(begin(fullEdges, Standard()), end(fullEdges, Standard()), LessMsaEdgeCargo_<TVertexDescriptor, TCargo>() );

    // Perform triplet extension
    typedef typename Iterator<TEdgeString, Standard>::Type TEdgeIter;
    TEdgeIter itEdges1 = begin(fullEdges, Standard());
    TEdgeIter itEdgesEnd = end(fullEdges, Standard());
    for(; itEdges1 != itEdgesEnd; ++itEdges1) {
        for(TEdgeIter itEdges2 = itEdges1; ++itEdges2 != itEdgesEnd;) {
            if ((*itEdges1).v1 != (*itEdges2).v1) break;
            if (sequenceId(g, (*itEdges1).v2) != sequenceId(g, (*itEdges2).v2)) {
                TCargo weight = ((*itEdges1).c < (*itEdges2).c) ? (*itEdges1).c : (*itEdges2).c;
                if ((*itEdges1).v2 < (*itEdges2).v2) {
                    TEdgeMapIter pos = newEMap.find(TNewEdge((*itEdges1).v2, (*itEdges2).v2));
                    if (pos != newEMap.end()) (*pos).second += weight;
                    else newEMap.insert(std::make_pair(TNewEdge((*itEdges1).v2, (*itEdges2).v2), weight));
                } else {
                    TEdgeMapIter pos = newEMap.find(TNewEdge((*itEdges2).v2, (*itEdges1).v2));
                    if (pos != newEMap.end()) (*pos).second += weight;
                    else newEMap.insert(std::make_pair(TNewEdge((*itEdges2).v2, (*itEdges1).v2), weight));
                }
            }
        }
    }
    clear(fullEdges);

    // Insert edges
    for(TEdgeMapIter itE = newEMap.begin(); itE != newEMap.end(); ++itE)
        addEdge(g, (*itE).first.first, (*itE).first.second, (*itE).second);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGuideTree, typename TSeqGroups, typename TGroupRoot, typename TSize>
inline void
_subTreeSearch(TGuideTree& guideTree,
               TSeqGroups& seqGroups,
               TGroupRoot& groupRoot,
               TSize minMembers)
{
    typedef typename Value<TSeqGroups>::Type TSeqGroup;
    typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;

    // Initialization
    TVertexDescriptor rootVertex = getRoot(guideTree);
    TSize nVertices = numVertices(guideTree);
    TSize nSeq = (nVertices / 2) + 1;

    // Number of subsequent leaves for each node
    typedef typename Iterator<TGuideTree, BfsIterator>::Type TBfsIterator;
    String<TSize> numLeaves;
    resizeVertexMap(numLeaves, guideTree);

    // All vertices in reversed bfs order
    typedef String<TVertexDescriptor> TVertexString;
    TVertexString vertices;
    resize(vertices, nVertices);

    // Walk through the tree in bfs order
    TBfsIterator bfsIt(guideTree, rootVertex);
    TSize pos = length(vertices) - 1;
    for(;!atEnd(bfsIt);goNext(bfsIt), --pos) {
        if (isLeaf(guideTree, *bfsIt)) property(numLeaves, *bfsIt) = 1;
        else property(numLeaves, *bfsIt) = 0;
        value(vertices, pos) = *bfsIt;
    }

    // Count the number of leaves for each internal node
    typedef typename Iterator<TVertexString, Standard>::Type TVertexIter;
    TVertexIter itVert = begin(vertices, Standard());
    TVertexIter itVertEnd = end(vertices, Standard());
    for(;itVert != itVertEnd; ++itVert) {
        if (!isLeaf(guideTree, *itVert)) {
            typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
            TAdjacencyIterator adjIt(guideTree, *itVert);
            for(;!atEnd(adjIt);goNext(adjIt)) property(numLeaves, *itVert) += property(numLeaves, *adjIt);
        }
    }

    // Delineate the groups
    itVert = begin(vertices, Standard());
    for(;itVert != itVertEnd; ++itVert) {
        if (property(numLeaves, *itVert) >= minMembers) {
            appendValue(seqGroups, TSeqGroup(), Generous());
            appendValue(groupRoot, *itVert, Generous());
            TSize elem = length(seqGroups) - 1;
            collectLeaves(guideTree, *itVert, seqGroups[elem]);
            property(numLeaves, *itVert) = 0;
            // Do not take any parent of the group root
            if (*itVert != rootVertex) {
                TVertexDescriptor pVert = parentVertex(guideTree, *itVert);
                while(pVert != rootVertex) {
                    property(numLeaves, pVert) = 0;
                    pVert = parentVertex(guideTree, pVert);
                }
                property(numLeaves, pVert) = 0;
            }
        }
    }
    if (!length(seqGroups)) {
        appendValue(seqGroups, TSeqGroup());
        appendValue(groupRoot, rootVertex);
        collectLeaves(guideTree, rootVertex, seqGroups[0]);
    }

    // Label all internal vertices with the closest root node
    typedef Pair<TSize, TSize> TDistGroup; // Distance, group index
    String<TDistGroup> closestRoot;
    resize(closestRoot, getIdUpperBound(_getVertexIdManager(guideTree)), TDistGroup(0,0), Exact());
    for(TSize i=0; i< (TSize) length(groupRoot); ++i) {
        TVertexDescriptor v = groupRoot[i];
        TSize dist = 0;
        while(v != rootVertex) {
            ++dist;
            v = parentVertex(guideTree, v);
            if ((property(closestRoot,v).i1 == 0) ||
                (property(closestRoot,v).i1 > dist)) {
                property(closestRoot, v) = TDistGroup(dist,i);
            }
        }
    }

    // Find ungrouped vertices
    typedef typename Iterator<TSeqGroup, Standard>::Type TSeqGroupIter;
    TSeqGroup allGroupedLeaves;
    for(TSize i=0; i< (TSize) length(seqGroups); ++i) {
        TSeqGroupIter itSeqGroup = begin(seqGroups[i], Standard());
        TSeqGroupIter itSeqGroupEnd = end(seqGroups[i], Standard());
        for(;itSeqGroup != itSeqGroupEnd; ++itSeqGroup)
            appendValue(allGroupedLeaves, *itSeqGroup, Generous());
    }
    appendValue(allGroupedLeaves, nSeq);
    std::sort(begin(allGroupedLeaves, Standard()), end(allGroupedLeaves, Standard()));
    TSeqGroupIter itSeqGroup = begin(allGroupedLeaves, Standard());
    TSeqGroupIter itSeqGroupEnd = end(allGroupedLeaves, Standard());
    TSize leftover = 0;
    TSeqGroup ungroupedLeaves;
    for(;itSeqGroup != itSeqGroupEnd; ++itSeqGroup, ++leftover) {
        while (leftover<*itSeqGroup) {
            appendValue(ungroupedLeaves, leftover, Generous());
            ++leftover;
        }
    }


    //std::cout << guideTree << std::endl;
    //std::cout << nSeq << std::endl;
    //std::cout << length(groupRoot) << std::endl;
    //std::cout << length(allGroupedLeaves) << std::endl;
    //std::cout << length(ungroupedLeaves) << std::endl;

    // Group the ungrouped vertices to the closest group
    clear(allGroupedLeaves);
    itSeqGroup = begin(ungroupedLeaves, Standard());
    itSeqGroupEnd = end(ungroupedLeaves, Standard());
    for(; itSeqGroup != itSeqGroupEnd; ++itSeqGroup) {
        TVertexDescriptor v = *itSeqGroup;
        while(v != rootVertex) {
            v = parentVertex(guideTree, v);
            if (property(closestRoot,v).i1 != 0) {
                appendValue(seqGroups[property(closestRoot,v).i2], *itSeqGroup, Generous());
                break;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TSize>
inline void
tripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                        TGuideTree& guideTree,
                        TSize minMembers)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename VertexDescriptor<TGuideTree>::Type TTreeVertex;
    TStringSet& strSet = stringSet(g);
    TSize nSeq = length(strSet);

    // Identify large subtrees
    typedef String<TSize> TSequenceGroups;
    String<TSequenceGroups> seqGroup;
    typedef String<TTreeVertex> TGroupRoot;
    TGroupRoot groupRoot;
    _subTreeSearch(guideTree, seqGroup, groupRoot, minMembers);

    // Label the subtree sequences
    String<TSize> seqLabels;
    resize(seqLabels, nSeq);
    typedef typename Iterator<TSequenceGroups, Standard>::Type TSeqSetIter;
    for(TSize i=0; i< (TSize) length(seqGroup); ++i) {
        TSeqSetIter itSeqGroup = begin(seqGroup[i], Standard());
        TSeqSetIter itSeqGroupEnd = end(seqGroup[i], Standard());
        for(;itSeqGroup != itSeqGroupEnd; ++itSeqGroup)
            seqLabels[*itSeqGroup] = i;
    }

    // Store all edges
    typedef std::pair<TVertexDescriptor, TVertexDescriptor> TNewEdge;
    typedef std::map<TNewEdge, TCargo> TEdgeMap;
    typedef typename TEdgeMap::iterator TEdgeMapIter;
    TEdgeMap newEMap;
    typedef MsaEdgeCargo_<TVertexDescriptor, TCargo> TEdge;
    typedef String<TEdge> TEdgeString;
    TEdgeString initialEdges;
    TEdgeIterator itE(g);
    for(;!atEnd(itE);goNext(itE)) {
        TVertexDescriptor sV = sourceVertex(itE);
        TVertexDescriptor tV = targetVertex(itE);
        TCargo c = cargo(*itE);
        appendValue(initialEdges, TEdge(sV, tV, c), Generous());
        newEMap.insert(std::make_pair(TNewEdge(sV, tV), c));
    }
    clearEdges(g);


    // Perform triplet extension
    typedef typename Iterator<TEdgeString, Standard>::Type TEdgeIter;
    TEdgeString fullEdges;
    for(TSize i=0; i< (TSize) length(seqGroup); ++i) {
        TEdgeIter itInitial = begin(initialEdges, Standard());
        TEdgeIter itInitialEnd = end(initialEdges, Standard());
        for(; itInitial != itInitialEnd; ++itInitial) {
            TSize seq1 = idToPosition(strSet, sequenceId(g, (*itInitial).v1));
            TSize seq2 = idToPosition(strSet, sequenceId(g, (*itInitial).v2));
            if ((seqLabels[seq1] == i) && (seqLabels[seq2] == i)) {
                appendValue(fullEdges, *itInitial, Generous());
                appendValue(fullEdges, TEdge((*itInitial).v2, (*itInitial).v1, (*itInitial).c), Generous());
            }
        }
        std::sort(begin(fullEdges, Standard()), end(fullEdges, Standard()), LessMsaEdgeCargo_<TVertexDescriptor, TCargo>() );
        typedef typename Iterator<TEdgeString, Standard>::Type TEdgeIter;
        TEdgeIter itEdges1 = begin(fullEdges, Standard());
        TEdgeIter itEdgesEnd = end(fullEdges, Standard());
        for(; itEdges1 != itEdgesEnd; ++itEdges1) {
            for(TEdgeIter itEdges2 = itEdges1; ++itEdges2 != itEdgesEnd;) {
                if ((*itEdges1).v1 != (*itEdges2).v1) break;
                if (sequenceId(g, (*itEdges1).v2) != sequenceId(g, (*itEdges2).v2)) {
                    TCargo weight = ((*itEdges1).c < (*itEdges2).c) ? (*itEdges1).c : (*itEdges2).c;
                    if ((*itEdges1).v2 < (*itEdges2).v2) {
                        TEdgeMapIter pos = newEMap.find(TNewEdge((*itEdges1).v2, (*itEdges2).v2));
                        if (pos != newEMap.end()) (*pos).second += weight;
                        else newEMap.insert(std::make_pair(TNewEdge((*itEdges1).v2, (*itEdges2).v2), weight));
                    } else {
                        TEdgeMapIter pos = newEMap.find(TNewEdge((*itEdges2).v2, (*itEdges1).v2));
                        if (pos != newEMap.end()) (*pos).second += weight;
                        else newEMap.insert(std::make_pair(TNewEdge((*itEdges2).v2, (*itEdges1).v2), weight));
                    }
                }
            }
        }
        clear(fullEdges);
    }

    // Insert edges
    for(TEdgeMapIter itE = newEMap.begin(); itE != newEMap.end(); ++itE)
        addEdge(g, (*itE).first.first, (*itE).first.second, (*itE).second);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
graphBasedTripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // Two tasks:
    // 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
    // 2) Augment all existing edges
    String<TCargo> newCargoMap;
    typedef String<TVertexDescriptor> TVertexString;
    typedef String<TCargo> TCargoString;
    TVertexString edges_vertices;
    TCargoString edges_cargo;

    // Triplet Extension
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // Remember the old cargo
    resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
    TEdgeIterator it(g);
    for(;!atEnd(it);++it)
        property(newCargoMap, *it) = cargo(*it);

    // Iterate over all vertices
    for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
        TOutEdgeIterator outIt1(g, *itVertex);
        while (!atEnd(outIt1)) {
            TOutEdgeIterator outIt2 = outIt1;
            goNext(outIt2);
            // Consider always 2 neighbors
            while (!atEnd(outIt2)) {
                TVertexDescriptor tV1 = targetVertex(outIt1);
                TVertexDescriptor tV2 = targetVertex(outIt2);
                if (sequenceId(g, tV1) != sequenceId(g,tV2)) {
                    TEdgeDescriptor e = findEdge(g, tV1, tV2);
                    if (e == 0) {
                        // New edge
                        TCargo val = (cargo(*outIt1) < cargo(*outIt2)) ? cargo(*outIt1) : cargo(*outIt2);

                        // Remember the edge with cargo
                        appendValue(edges_vertices, tV1, Generous());
                        appendValue(edges_vertices, tV2, Generous());
                        appendValue(edges_cargo, val, Generous());
                    } else {
                        // Increase weight of existing edge
                        if (cargo(*outIt2) > cargo(*outIt1)) property(newCargoMap, e) += cargo(*outIt1);
                        else property(newCargoMap, e) += cargo(*outIt2);
                    }
                }
                goNext(outIt2);
            }
            goNext(outIt1);
        }
    }


    // Assign the new weights and clean-up the cargo map
    TEdgeIterator itE(g);
    for(;!atEnd(itE);++itE) cargo(*itE) = property(newCargoMap, *itE);
    clear(newCargoMap);

    // Add edges
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TVertexString, Standard>::Type TVertexStringIter;
    typedef typename Iterator<TCargoString, Standard>::Type TCargoStringIter;

    // Finally add the new edges created by the triplet approach
    TVertexStringIter itV = begin(edges_vertices, Standard());
    TVertexStringIter itVNext = begin(edges_vertices, Standard());
    ++itVNext;
    TVertexStringIter endIt = end(edges_vertices, Standard());
    TCargoStringIter itC = begin(edges_cargo, Standard());
    for(;itV != endIt; itV += 2, itVNext +=2, ++itC) {
        // The same edge could have been created multiple times, so check if it exists
        TEdgeDescriptor e = findEdge(g, *itV, *itVNext);
        if (e == 0) addEdge(g, *itV, *itVNext, *itC);
        else cargo(e) += *itC;
    }
    SEQAN_ASSERT(itC == end(edges_cargo));
}

////////////////////////////////////////////////////////////////////////////////
//
//template<typename TStringSet, typename TCargo, typename TSpec>
//inline void
//reducedTripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
//{
//    SEQAN_CHECKPOINT
//    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
//    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
//    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
//    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
//    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
//    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
//
//
//    // Just augment existing edges
//    String<TCargo> newCargoMap;
//    resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
//    TEdgeIterator it(g);
//    for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));
//
//    // Iterate over all vertices
//    for(TVertexIterator itVertex(g);!atEnd(itVertex);++itVertex) {
//        TOutEdgeIterator outIt1(g, *itVertex);
//        while (!atEnd(outIt1)) {
//            TOutEdgeIterator outIt2 = outIt1;
//            goNext(outIt2);
//            // Consider always 2 neighbors
//            while (!atEnd(outIt2)) {
//                TVertexDescriptor tV1 = targetVertex(outIt1);
//                TVertexDescriptor tV2 = targetVertex(outIt2);
//                if (sequenceId(g, tV1) != sequenceId(g,tV2)) {
//                    TEdgeDescriptor e = findEdge(g, tV1, tV2);
//                    if (e != 0) {
//                        // Increase weight of existing edge
//                        if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
//                        else property(newCargoMap, e) += getCargo(*outIt2);
//                    }
//                }
//                goNext(outIt2);
//            }
//            goNext(outIt1);
//        }
//    }
//
//    // Assign the new weights and clean-up the cargo map
//    TEdgeIterator itE(g);
//    for(;!atEnd(itE);goNext(itE)) cargo(value(itE)) = property(newCargoMap, value(itE));
//    clear(newCargoMap);
//}

//////////////////////////////////////////////////////////////////////////////
// Sum of Pairs Scoring
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// This version is sensitive to gap openings

/*!
 * @fn sumOfPairsScore
 * @headerfile <seqan/graph_msa.h>
 * @brief Given a multiple alignment, this function calculates the sum-of-pairs score.
 *
 * @signature TScoreValue sumOfPairsScore(graph, scoringScheme);
 *
 * @param[in] graph         An @link AlignmentGraph @endlink to use for the evaluation.
 * @param[in] ScoringScheme The @link Score @endlink to use.
 *
 * @return TScoreValue The SOP score of the MSA  (Metafunction: @link Score#Value @endlink of the type of
 *                     <tt>scoringScheme</tt>).
 *
 * This function does NOT assume independent columns.  That is, gap openings are properly scored.  If you want the fast
 * version assuming independ columns use sumOfPairsScoreInd.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline typename Value<TScore>::Type
sumOfPairsScore(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                TScore const& score_type)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;

    // Convert the graph
    String<char> mat;
    convertAlignment(g, mat);
    char gapChar = gapValue<char>();

    TScoreValue gap = scoreGapExtend(score_type);
    TScoreValue gapOpen = scoreGapOpen(score_type);
    TSize nseq = length(stringSet(g));
    TSize len = length(mat) / nseq;

    bool gapOpeni = false;
    bool gapOpenj = false;
    TScoreValue totalScore = 0;
    for(TSize i = 0; i<nseq-1; ++i) {
        for(TSize j=i+1; j<nseq; ++j) {
            for(TSize k=0;k<len; ++k) {
                if (value(mat, i*len+k) != gapChar) {
                    if (value(mat, j*len + k) != gapChar) {
                        gapOpeni = false;
                        gapOpenj = false;
                        totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
                    } else {
                        if (gapOpenj) {
                            totalScore += gap;
                        } else {
                            gapOpenj = true;
                            totalScore += gapOpen;
                        }
                    }
                } else if (value(mat, j*len + k) != gapChar) {
                        if (gapOpeni) {
                            totalScore += gap;
                        } else {
                            gapOpeni = true;
                            totalScore += gapOpen;
                        }
                }
            }
        }
    }
    return totalScore;
}

// TODO(holtgrew): Document sumOfPairsScoreInd.

//////////////////////////////////////////////////////////////////////////////
// This version is insensitive to gap openings, assumes independent columns
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline typename Value<TScore>::Type
sumOfPairsScoreInd(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                   TScore const& score_type)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;

    // Convert the graph
    String<char> mat;
    convertAlignment(g, mat);
    char gapChar = gapValue<char>();

    TScoreValue gap = scoreGapExtend(score_type);
    TSize nseq = length(stringSet(g));
    TSize len = length(mat) / nseq;

    TScoreValue totalScore = 0;
    for(TSize k=0;k<len; ++k) {
        for(TSize i = 0; i<nseq-1; ++i) {
            for(TSize j=i+1; j<nseq; ++j) {
                if (value(mat, i*len+k) != gapChar) {
                    if (value(mat, j*len + k) != gapChar) {
                        totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
                    } else totalScore += gap;
                } else if (value(mat, j*len + k) != gapChar) {
                        totalScore += gap;
                }
            }
        }
    }
    return totalScore;
}


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn alignmentEvaluation
 * @headerfile <seqan/graph_msa.h>
 * @brief Given a multiple sequence alignment, this function calculates different kinds of alignment statistics.
 *
 * @signature TScoreVal alignmentEvaluation(graph, scoringScheme, gapExCount, gapCount, pairCount, numPairs, len);
 *
 * @param[in]  graph         The @link AlignmentGraph @endlink object to compute statistics for.
 * @param[in]  scoringScheme The @link Score @endlink object to use.
 * @param[out] gapExCount    The number of gap extensions.
 * @param[out] gapCount      The number of gaps.
 * @param[out] pairCount     The number of aligned pairs.
 * @param[out] numPairs      Counter for each pair. A @link String @endlink resized to size <tt>n*n</tt>, where
                             <tt>n</tt> is the @link FiniteOrderedAlphabetConcept#ValueSize @endlink of the alphabet
                             of the aligned sequences.
 * @param[out] len           Alignment length.
 *
 * @return TScoreVal The score of the alignment  (Metafunction: @link Score#Value @endlink).
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TSize>
inline typename Value<TScore>::Type
alignmentEvaluation(Graph<Alignment<TStringSet, TCargo, TSpec> > const & g,
                    TScore const& score_type,
                    TSize& gapExCount,
                    TSize& gapCount,
                    TSize& pairCount,
                    String<TSize> & numPairs,
                    TSize& len)
{
    //typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;

    // Initialization;
    gapExCount = 0;
    gapCount = 0;
    pairCount = 0;
    clear(numPairs);

    // Convert the graph
    String<char> mat;
    convertAlignment(g, mat);
    char gapChar = gapValue<char>();

    TScoreValue gap = scoreGapExtend(score_type);
    TScoreValue gapOpen = scoreGapOpen(score_type);
    TSize nseq = length(stringSet(g));
    len = length(mat) / nseq;

    bool gapOpeni = false;
    bool gapOpenj = false;
    TScoreValue totalScore = 0;
    resize(numPairs, alphSize * alphSize, 0);
    for(TSize i = 0; i<nseq-1; ++i) {
        for(TSize j=i+1; j<nseq; ++j) {
            for(TSize k=0;k<len; ++k) {
                if (value(mat, i*len+k) != gapChar) {
                    if (value(mat, j*len + k) != gapChar) {
                        gapOpeni = false;
                        gapOpenj = false;
                        ++pairCount;
                        TSize index1 = ordValue(TAlphabet(value(mat, i*len+k)));
                        TSize index2 = ordValue(TAlphabet(value(mat, j*len + k)));
                        value(numPairs, index1 * alphSize + index2) += 1;
                        totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
                    } else {
                        if (gapOpenj) {
                            ++gapExCount;
                            totalScore += gap;
                        } else {
                            gapOpenj = true;
                            ++gapCount;
                            totalScore += gapOpen;
                        }
                    }
                } else if (value(mat, j*len + k) != gapChar) {
                        if (gapOpeni) {
                            ++gapExCount;
                            totalScore += gap;
                        } else {
                            ++gapCount;
                            gapOpeni = true;
                            totalScore += gapOpen;
                        }
                }
            }
        }
    }
    return totalScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSource, typename TSpec2>
inline bool
convertAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& gAlign,
                 Align<TSource, TSpec2>& align)
{
    // Pipe into Align data structure
    String<char> mat;
    if (convertAlignment(gAlign, mat)) {
        typedef Align<TSource, TSpec2> TAlign;
        typedef typename Size<TAlign>::Type TSize;
        typedef typename Row<TAlign>::Type TRow;
        typedef typename Iterator<TRow>::Type TRowIterator;
        TStringSet& sourceSet = stringSet(gAlign);
        TSize nseq = length(sourceSet);
        clearGaps(align);
        if (empty(rows(align))) {
            resize(rows(align), nseq);
            for(TSize i = 0; i<nseq; ++i) assignSource(row(align, i), sourceSet[i]);
        }
        String<TRowIterator> rowIter;
        resize(rowIter, nseq);
        for(TSize i = 0; i<nseq; ++i) value(rowIter, i) = begin(row(align, i));
        TSize lenMat = length(mat);
        TSize colLen = lenMat / nseq;
        TSize gapCount = 0;
        char gapChar = gapValue<char>();
        for(TSize alignRow = 0; alignRow < nseq; ++alignRow) {
            for(TSize pos = alignRow * colLen; pos < (alignRow + 1) * colLen; ++pos) {
                if (value(mat, pos) != gapChar) {
                    if (gapCount) {
                        insertGaps(value(rowIter, alignRow), gapCount);
                        goFurther(value(rowIter, alignRow), gapCount);
                        gapCount = 0;
                    }
                    goNext(value(rowIter,alignRow));
                } else ++gapCount;
            }
            if (gapCount) {
                insertGaps(value(rowIter, alignRow), gapCount);
                goFurther(value(rowIter, alignRow), gapCount);
                gapCount = 0;
            }
        }
    } else return false;
    return true;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSource, typename TSpec2, typename TStringSet, typename TCargo, typename TSpec>
inline bool
convertAlignment(Align<TSource, TSpec2> const& align,
                 Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign)
{
    typedef Align<TSource, TSpec2> const TAlign;
    //typedef typename Value<TSource>::Type TAlphabet;
    typedef typename Size<TAlign>::Type TSize;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;
    clearVertices(gAlign);
    if (!length(stringSet(gAlign))) {
        TStringSet sourceSet = stringSet(const_cast<Align<TSource, TSpec2>&>(align));
        assignStringSet(gAlign, sourceSet);
    }
    TSize nseq = length(rows(align));
    String<TRowIterator> rowIter;
    String<TRowIterator> rowIterEnd;
    resize(rowIter, nseq);
    resize(rowIterEnd, nseq);
    for(TSize i = 0; i<nseq; ++i) {
        value(rowIter, i) = begin(row(align, i), Standard());
        value(rowIterEnd, i) = end(row(align, i), Standard());
    }
    String<Fragment<> > matches;
    for(TSize alignRow1 = 0; alignRow1 < nseq; ++alignRow1) {
        for(TSize alignRow2 = alignRow1 + 1; alignRow2 < nseq; ++alignRow2) {
            TRowIterator pos1 = value(rowIter,alignRow1);
            TRowIterator pos2 = value(rowIter,alignRow2);
            TSize alignPos = 0;
            TSize length = 0;
            TSize offset1 = 0;
            TSize offset2 = 0;
            for(;pos1 != value(rowIterEnd, alignRow1); ++pos1, ++pos2, ++alignPos) {
                if ((isGap(pos1)) || (isGap(pos2))) {
                    if (length) {
                        appendValue(matches, Fragment<>(alignRow1, alignPos - offset1 - length, alignRow2, alignPos - offset2 - length, length));
                        length = 0;
                    }
                    if (isGap(pos1)) ++offset1;
                    if (isGap(pos2)) ++offset2;
                } else ++length;
            }
            if (length) appendValue(matches, Fragment<>(alignRow1, alignPos - offset1 - length, alignRow2, alignPos - offset2 - length, length));
        }
    }
    matchRefinement(matches,stringSet(gAlign),gAlign);
    return true;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
