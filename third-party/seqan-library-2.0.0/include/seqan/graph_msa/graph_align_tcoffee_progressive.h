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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_PROGRESSIVE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Progressive Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPosition, typename TSequence>
inline void
_buildLeafString(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                 TPosition const pos,
                 TSequence& alignSeq)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Id<TGraph>::Type TId;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TSequence>::Type TVertexString;

    //TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    TStringSet& str = stringSet(g);
    TId seqId = positionToId(str, pos);
    TSize lenRoot = length(str[pos]);
    TSize i = 0;
    while(i<lenRoot) {
        TVertexDescriptor nextVertex = findVertex(const_cast<TGraph&>(g), seqId, i);
        //SEQAN_ASSERT(nextVertex != nilVertex);
        //if (nextVertex == nilVertex) {
        //    std::cout << "Warning: Nil Vertex" << std::endl;
        //    TSize j = i + 1;
        //    while ((j < lenRoot) && (findVertex(const_cast<TGraph&>(g), seqId, j) == nilVertex)) ++j;
        //    nextVertex = addVertex(const_cast<TGraph&>(g), seqId, i, j-i);
        //}
        TVertexString vs;
        appendValue(vs, nextVertex);
        appendValue(alignSeq, vs, Generous());
        i += fragmentLength(g, nextVertex);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSegmentString, typename TOutGraph>
inline void
_createAlignmentGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                      TSegmentString& alignSeq,
                      TOutGraph& gOut)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<TVertexDescriptor> TVertexString;

    // Create the alignment graph
    TSize alignSeqLen = length(alignSeq);
    for(TSize i = 0; i<alignSeqLen;++i) {
        TVertexString& alignSeq_i = alignSeq[i];
        TSize len_i = length(alignSeq_i);
        for(TSize j=0; j<len_i; ++j) {
            TVertexDescriptor v = alignSeq_i[j];
            SEQAN_ASSERT(fragmentBegin(g,v) < length(getValueById(stringSet(g), sequenceId(g,v))));
            SEQAN_ASSERT(fragmentLength(g,v) > 0);
            SEQAN_ASSERT(fragmentBegin(g,v) + fragmentLength(g,v) <= length(getValueById(stringSet(g), sequenceId(g,v))));
            TVertexDescriptor l = addVertex(gOut, sequenceId(g, v), fragmentBegin(g,v), fragmentLength(g,v));
            //std::cout << l << label(gOut, l) << ',';
            TSize count = 1;
            for(TSize k = j; k>0; --k) {
                //SEQAN_ASSERT(fragmentLength(gOut,l) == fragmentLength(gOut,l - count));
                addEdge(gOut, (TVertexDescriptor) (l - count), (TVertexDescriptor) l);
                ++count;
            }
        }
        //std::cout << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn progressiveAlignment
 * @headerfile <seqan/graph_msa.h>
 * @brief Perform a progressive multiple sequence alignment (MSA).
 *
 * @signature void progressiveAlignment(inputGraph, guideTree, outputGraph);
 *
 * @param[in]  inputGraph  A @link AlignmentGraph @endlink with multiple sequence information.
 * @param[in]  guideTree   A @link Tree @endlink to use as the guide tree.
 * @param[out] outputGraph An @link AlignmentGraph @endlink for the final MSA.
 */

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph>
inline void
progressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                     TGuideTree& tree,
                     TOutGraph& gOut)
{
    SEQAN_CHECKPOINT
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
    typedef typename Iterator<TGuideTree, BfsIterator>::Type TBfsIterator;
    typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
    typedef String<TVertexDescriptor> TVertexString;
    typedef String<TVertexString> TSegmentString;

    // Initialization
    TVertexDescriptor rootVertex = getRoot(tree);
    TSize nVertices = numVertices(tree);

    // Vertices in reversed bfs order
    TVertexString vertices;
    resize(vertices, nVertices);

    // All Strings of Strings of vertices for each node of the guide tree
    String<TSegmentString> segString;
    resize(segString, nVertices);

    // Walk through the tree in bfs order
    typedef typename Iterator<TVertexString, Standard>::Type TVertexIter;
    TVertexIter itVert = begin(vertices, Standard());
    TVertexIter itVertEnd = end(vertices, Standard());
    --itVertEnd;
    TBfsIterator bfsIt(tree, rootVertex);
    for(;!atEnd(bfsIt);goNext(bfsIt), --itVertEnd)
        *itVertEnd = *bfsIt;

    // Progressive alignment
    itVert = begin(vertices, Standard());
    itVertEnd = end(vertices, Standard());
    for(;itVert != itVertEnd; ++itVert) {
        if(isLeaf(tree, *itVert)) _buildLeafString(g, *itVert, segString[*itVert]);
        else {
            // Align the two children (Binary tree)
            TAdjacencyIterator adjIt(tree, *itVert);
            TVertexDescriptor child1 = *adjIt; goNext(adjIt);
            heaviestCommonSubsequence(g, segString[child1], segString[*adjIt], segString[*itVert]);
            clear(segString[child1]);
            clear(segString[*adjIt]);
        }
    }

    // Create the alignment graph
    _createAlignmentGraph(g, segString[rootVertex], gOut);
}



//////////////////////////////////////////////////////////////////////////////
// Progressive Matching, just for testing purposes
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSegmentString, typename TEdgeMap, typename TOutGraph>
inline void
_createMatchingGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                     TSegmentString& alignSeq,
                     TEdgeMap& edgeMap,
                     TOutGraph& gOut,
                     TEdgeMap& edgeMapOut)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef String<TVertexDescriptor> TVertexString;

    // Initialization
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

    // Create the matching graph
    clear(edgeMapOut);
    TSize alignSeqLen = length(alignSeq);
    for(TSize i = 0; i<alignSeqLen;++i) {
        TVertexString& alignSeq_i = alignSeq[i];
        TSize len_i = length(alignSeq_i);
        for(TSize j=0; j<len_i - 1; ++j) {
            for(TSize k=j+1; k<len_i; ++k) {
                TVertexDescriptor v1 = getValue(alignSeq_i, j);
                TVertexDescriptor v2 = getValue(alignSeq_i, k);
                if ((v1 == nilVertex) || (v2 == nilVertex)) continue;

                TVertexDescriptor v1New = findVertex(gOut, sequenceId(g, v1), fragmentBegin(g,v1));
                if (v1New == nilVertex) v1New = addVertex(gOut, sequenceId(g, v1), fragmentBegin(g,v1), fragmentLength(g,v1));
                TVertexDescriptor v2New = findVertex(gOut, sequenceId(g, v2), fragmentBegin(g,v2));
                if (v2New == nilVertex) v2New = addVertex(gOut, sequenceId(g, v2), fragmentBegin(g,v2), fragmentLength(g,v2));

                TEdgeDescriptor e = findEdge(g, v1, v2);
                addEdge(gOut, v1New, v2New, cargo(e));
                appendValue(edgeMapOut, property(edgeMap, e));
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString, typename TOutString>
inline TCargo
heaviestMatching(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
                 TString const& str1,
                 TString const& str2,
                 TOutString& align)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef __int64 TLargeSize;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TString>::Type TVertexSet;

    TSize m = length(str1);  // How many sets of vertex descriptors in seq1
    TSize n = length(str2);  // How many sets of vertex descriptors in seq1

    // Size of the sequences
    // Note for profile alignments every member of the sequence is a String!!! of vertex descriptors
    TCargo divider = (TCargo) length(str1[0]) * (TCargo) length(str2[0]);
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

    // Fill the vertex to position map for str1
    // Remember for each vertex descriptor the position in the sequence
    typedef std::map<TVertexDescriptor, TSize> TVertexToPosMap;
    typedef typename TVertexToPosMap::const_iterator TVertexToPosMapIter;
    TVertexToPosMap map;
    typedef typename Iterator<TString const>::Type TStringIterConst;
    typedef typename Iterator<TVertexSet const>::Type TVertexSetIterConst;
    TStringIterConst itStrEnd1 = end(str1);
    TSize pos = 0;
    for(TStringIterConst itStr1 = begin(str1);itStr1 != itStrEnd1;++itStr1, ++pos) {
        TVertexSetIterConst itVEnd = end(getValue(itStr1));
        for(TVertexSetIterConst itV = begin(getValue(itStr1));itV != itVEnd;++itV) {
            if (*itV != nilVertex) map.insert(std::make_pair(*itV, pos));
        }
    }

    // Given a full graph, what positions are occupied?
    typedef std::set<TLargeSize> TOccupiedPositions;
    TOccupiedPositions occupiedPositions;
    TStringIterConst itStrEnd2 = end(str2);
    TSize posItStr2 = 0;
    for(TStringIterConst itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2, ++posItStr2) {
        TVertexSetIterConst itVEnd = end(getValue(itStr2));
        for(TVertexSetIterConst itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
            if (*itV != nilVertex) {
                TOutEdgeIterator itOut(g, *itV);
                for(;!atEnd(itOut); ++itOut) {
                    // Target vertex must be in the map
                    TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
                    if (pPos != map.end()) occupiedPositions.insert( (TLargeSize) (pPos->second) * (TLargeSize) n + (TLargeSize) (n - posItStr2 - 1) );
                }
            }
        }
    }
    // Map the occupied positions to slots
    typedef std::map<TLargeSize, TSize> TPositionToSlotMap;
    TPositionToSlotMap posToSlotMap;
    TSize counter = 0;
    for(typename TOccupiedPositions::const_iterator setIt = occupiedPositions.begin();setIt != occupiedPositions.end(); ++setIt, ++counter) {
        posToSlotMap.insert(std::make_pair(*setIt, counter));
    }
    occupiedPositions.clear();

    // Walk through str2 and fill in the weights of the actual edges
    typedef String<TCargo> TWeights;
    typedef typename Iterator<TWeights>::Type TWeightsIter;
    TWeights weights;
    resize(weights, posToSlotMap.size(),0);
    posItStr2 = 0;
    for(TStringIterConst itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2, ++posItStr2) {
        TVertexSetIterConst itVEnd = end(getValue(itStr2));
        for(TVertexSetIterConst itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
            if (*itV != nilVertex) {
                TOutEdgeIterator itOut(g, *itV);
                for(;!atEnd(itOut); ++itOut) {
                    // Target vertex must be in the map
                    TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
                    if (pPos != map.end()) weights[posToSlotMap[ (TLargeSize) (pPos->second) * (TLargeSize) n + (TLargeSize) (n - posItStr2 - 1) ]] += (TCargo) cargo(*itOut);
                }
            }
        }
    }
    map.clear();
    // Average weights
    TWeightsIter itWeights = begin(weights);
    TWeightsIter itWeightsEnd = begin(weights);
    for(;itWeights != itWeightsEnd; ++itWeights) *itWeights /= divider;

    // Create the match graph
    Graph<Undirected<void> > matchGraph;
    typedef typename VertexDescriptor<Graph<Undirected<void> > >::Type TVD;
    String<bool> vertexMap;
    resize(vertexMap, m + n, false);
    for(TSize i = 0; i < m; ++i) addVertex(matchGraph);
    for(TSize i = 0; i < n; ++i) value(vertexMap, addVertex(matchGraph)) = true;
    for(typename TPositionToSlotMap::const_iterator mapIt = posToSlotMap.begin();mapIt != posToSlotMap.end(); ++mapIt) {
        TLargeSize pos = mapIt->first;
        TSize i = (TSize) (pos / (TLargeSize) n);   // Get the index in str1
        TSize j = n - 1 - (TSize) (pos % (TLargeSize) n); // Get the index in str2
        addEdge(matchGraph, i, m+j);
    }

    // Compute the best matching
    typedef String<Pair<TVD, TVD> > TEdges;
    TEdges edges;

    //std::fstream strm;
    //strm.open("Z:\\my_graph.dot", std::ios_base::out | std::ios_base::trunc);
    //writeRecords(strm,matchGraph,DotDrawing());
    //strm.close();

    TCargo val = weightedBipartiteMatching(matchGraph, vertexMap, weights, edges);

    // Retrieve the aligned segments
    TSize seqsInStr1 = length(str1[0]);
    TSize seqsInStr2 = length(str2[0]);
    typedef typename Value<TString>::Type TVertexSet;
    //typedef typename Iterator<TString const, Rooted>::Type TStringIter;
    //typedef typename Iterator<TString, Rooted>::Type TSIter;
    typedef typename Iterator<TVertexSet const, Rooted>::Type TVertexSetIter;
    //typedef typename Iterator<TVertexSet, Rooted>::Type TIter;
    clear(align);
    // Retrieve all matches
    String<bool> matchedVertices;
    resize(matchedVertices, n + m, false);
    typedef typename Iterator<TEdges>::Type TEdgesIter;
    TEdgesIter itEdges = begin(edges);
    TEdgesIter itEdgesEnd = end(edges);
    for(;itEdges != itEdgesEnd; ++itEdges) {
        TSize i = (value(itEdges)).i1;
        TSize j = (value(itEdges)).i2 - m;
        value(matchedVertices, i) = true;
        value(matchedVertices, m+j) = true;
        TVertexSet tmp;
        resize(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
        TVertexSetIter itVEnd = end(value(str1,i));
        TSize count = 0;
        for(TVertexSetIter itV = begin(value(str1,i));itV != itVEnd;++itV) {
            tmp[count] = *itV;
            ++count;
        }
        TVertexSetIter itVEnd2 = end(value(str2,j));
        for(TVertexSetIter itV2 = begin(value(str2,j));itV2 != itVEnd2;++itV2) {
            tmp[count] = *itV2;
            ++count;
        }
        appendValue(align, tmp);
    }
    typedef typename Iterator<Graph<Undirected<void> >, VertexIterator>::Type TVertexIterator;
    TVertexIterator itVertex(matchGraph);
    for(;!atEnd(itVertex);++itVertex) {
        if (getProperty(matchedVertices, *itVertex) == true) continue;
        if (*itVertex < m) {
            TSize i = *itVertex;
            TVertexSet tmp;
            resize(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
            TVertexSetIter itVEnd = end(value(str1, i));
            TSize count = 0;
            for(TVertexSetIter itV = begin(value(str1, i));itV != itVEnd;++itV) {
                tmp[count] = *itV;
                ++count;
            }
            appendValue(align, tmp);
        } else {
            TSize j = *itVertex - m;
            TVertexSet tmp;
            resize(tmp, (seqsInStr1 + seqsInStr2), nilVertex);
            TVertexSetIter itVEnd = end(value(str2, j));
            TSize count = 0;
            for(TVertexSetIter itV = begin(value(str2, j));itV != itVEnd;++itV) {
                tmp[count] = *itV;
                ++count;
            }
            appendValue(align, tmp);
        }
    }

    return val;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence>
inline void
_recursiveProgressiveMatching(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                              TGuideTree& tree,
                              TVertexDescriptor const root,
                              TSequence& alignSeq)
{
    //typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    //typedef typename Size<TGraph>::Type TSize;
    //typedef typename Id<TGraph>::Type TId;
    typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
    //typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

    if(isLeaf(tree, root)) {
        _buildLeafString(g, root, alignSeq);
    } else {
        // Align the two children (Binary tree)
        typedef String<String<TVertexDescriptor> > TSegmentString;
        TSegmentString seq1;
        TSegmentString seq2;
        TAdjacencyIterator adjIt(tree, root);
        _recursiveProgressiveMatching(g,tree, *adjIt, seq1);
        goNext(adjIt);
        _recursiveProgressiveMatching(g,tree, *adjIt, seq2);

        heaviestMatching(g, seq1, seq2, alignSeq);

        //// Debug Code
        //for(TSize i = 0; i<length(alignSeq);++i) {
        //    std::cout << '(';
        //    for(TSize j=0; j<length(alignSeq[i]);++j) {
        //        std::cout << getValue(alignSeq[i], j) << ',';
        //    }
        //    std::cout << ')' << ',';
        //}
        //std::cout << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TEdgeMap, typename TOutGraph>
inline void
progressiveMatching(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
                    TGuideTree& tree,
                    TEdgeMap& edgeMap,
                    TOutGraph& gOut,
                    TEdgeMap& edgeMapOut)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    //typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<TVertexDescriptor> TVertexString;
    typedef String<TVertexString> TSegmentString;

    // Perform progressive alignment
    TSegmentString alignSeq;
    _recursiveProgressiveMatching(g,tree,getRoot(tree),alignSeq);

    // Create the alignment graph
    _createMatchingGraph(g, alignSeq, edgeMap, gOut, edgeMapOut);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
