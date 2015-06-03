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
// Author: Tobias Raussch <rausch@embl.de>
// ==========================================================================
// Implementation of Weighted-Bipartite-Matching algorithm.
//
// WARNING: Functionality not carefully tested!
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_WEIGHTED_BIPARTITE_MATCHING_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_WEIGHTED_BIPARTITE_MATCHING_H_

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
// Function weightedBipartiteMatching()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TVertexMap, typename TWeightMap, typename TEdges>
inline typename Value<TWeightMap>::Type
_weightedBipartiteMatching(Graph<TSpec>& g,
                              TVertexMap& vertMap,
                              TWeightMap& weightMap,
                              String<TEdges>& edges)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIter;
    typedef typename Value<TWeightMap>::Type TCargo;

    TSize numVert = numVertices(g);
    TCargo maxEdgeVal = 0;

    // Find an initial labeling
    String<TCargo> label;
    resizeVertexMap(label, g);
    TVertexIterator itV(g);
    for(;!atEnd(itV); goNext(itV)) {
        if (getProperty(vertMap, value(itV)) == true) value(label, value(itV)) = 0;
        else {
            TCargo maxCargo = 0;
            for(TOutEdgeIter itOutE(g, value(itV));!atEnd(itOutE); goNext(itOutE)) {
                if (property(weightMap, (value(itOutE))) > maxCargo) maxCargo = property(weightMap, (value(itOutE)));
            }
            value(label, value(itV)) = maxCargo;
            if (maxCargo > maxEdgeVal) maxEdgeVal = maxCargo;
        }
    }

    // Generate Equality Graph
    typedef Graph<Directed<void> > TEqualityGraph;
    typedef typename EdgeType<TEqualityGraph>::Type TEdgeStump;
    TEqualityGraph equalGraph;
    resize(equalGraph.data_vertex, length(_getVertexString(g)), (TEdgeStump*) 0);
    equalGraph.data_id_managerV = g.data_id_managerV;
    TEdgeIterator itE(g);
    for(;!atEnd(itE); goNext(itE)) {
        if (property(weightMap, (value(itE))) == property(label, sourceVertex(itE)) + property(label, targetVertex(itE))) {
            // For the Ford-Fulkerson all edges must go from true to false
            if (getProperty(vertMap, sourceVertex(itE)) == true) addEdge(equalGraph, targetVertex(itE), sourceVertex(itE));
            else addEdge(equalGraph, sourceVertex(itE), targetVertex(itE));
        }
    }

    // Find an initial bipartite matching
    clear(edges);
    TSize matchSize = bipartiteMatching(equalGraph, vertMap, edges);


    String<bool> free;
    String<TVertexDescriptor> reverseMatchMap;
    typedef std::set<TVertexDescriptor> TVertexSet;
    TVertexSet setS;
    TVertexSet setNeighborS;
    TVertexSet setT;
    while (matchSize != numVert / 2) {

        // Initialization
        setS.clear();
        setT.clear();
        setNeighborS.clear();
        clear(free);
        resize(free, getIdUpperBound(_getVertexIdManager(g)), true);
        clear(reverseMatchMap);
        resizeVertexMap(reverseMatchMap, g);

        // Find free vertex
        typedef typename Iterator<String<TEdges> >::Type TStringEdgeIter;
        TStringEdgeIter itSE = begin(edges);
        TStringEdgeIter itSEEnd = end(edges);
        for(;itSE != itSEEnd; goNext(itSE)) {
            value(free, (value(itSE)).i1) = false;
            value(free, (value(itSE)).i2) = false;
            value(reverseMatchMap, (value(itSE)).i2) = (value(itSE)).i1;
        }
        TVertexIterator itVert(g);
        for(;!atEnd(itVert); goNext(itVert)) {
            if ((getProperty(vertMap, value(itVert)) == false) &&
                (value(free, value(itVert)) == true)) {
                    setS.insert(value(itVert));
                    typedef typename Iterator<TEqualityGraph, OutEdgeIterator>::Type TOutEdgeIterator;
                    TOutEdgeIterator itOE(equalGraph, value(itVert));
                    for(;!atEnd(itOE); ++itOE) {
                        setNeighborS.insert(targetVertex(itOE));
                        if (value(free, targetVertex(itOE)) == true) setT.insert(targetVertex(itOE));
                    }
                    break;
            }
        }

        // Find matched vertices
        typedef typename TVertexSet::iterator TVertexSetIter;
        while (setNeighborS != setT) {
            TVertexSet diffSet;
            TVertexSetIter itT = setT.begin();
            TVertexSetIter itTEnd = setT.end();
            TVertexSetIter itN = setNeighborS.begin();
            TVertexSetIter itNEnd = setNeighborS.end();
            while (itN != itNEnd) {
                if ((itT == itTEnd) || (*itN < *itT)) { diffSet.insert(*itN); ++itN; }
                else { ++itN; ++itT; }
            }
            TVertexDescriptor y = *(diffSet.begin());
            setT.insert(y);
            setS.insert(value(reverseMatchMap, y));
            typedef typename Iterator<TEqualityGraph, OutEdgeIterator>::Type TOutEdgeIterator;
            TOutEdgeIterator itOE(equalGraph, value(reverseMatchMap, y));
            for(;!atEnd(itOE); ++itOE) {
                setNeighborS.insert(targetVertex(itOE));
                if (value(free, targetVertex(itOE)) == true) setT.insert(targetVertex(itOE));
            }
        }
        clear(reverseMatchMap);

        // Update Labels
        TCargo minVal = maxEdgeVal;
        TEdgeIterator itEdge(g);
        for(;!atEnd(itEdge); goNext(itEdge)) {
            TVertexDescriptor sV = sourceVertex(itEdge);
            TVertexDescriptor tV = targetVertex(itEdge);
            if (property(vertMap, sV) == true) {    TVertexDescriptor tmp = sV; sV = tV; tV = tmp;  }
            if ((setS.find(sV) != setS.end()) &&
                (setT.find(tV) == setT.end())) {
                TCargo thisVal = getProperty(label, sV) + getProperty(label, tV) - getProperty(weightMap, (value(itEdge)));
                if (thisVal < minVal) minVal = thisVal;
            }
        }
        TVertexIterator myVertexIt(g);
        for(;!atEnd(myVertexIt); goNext(myVertexIt)) {
            if (setS.find(value(myVertexIt)) != setS.end()) value(label, value(myVertexIt)) -= minVal;
            else if (setT.find(value(myVertexIt)) != setT.end()) value(label, value(myVertexIt)) += minVal;
        }

        // Build new equal graph
        clear(equalGraph);
        resize(equalGraph.data_vertex, length(_getVertexString(g)), (TEdgeStump*) 0);
        equalGraph.data_id_managerV = g.data_id_managerV;
        TEdgeIterator itE(g);
        for(;!atEnd(itE); goNext(itE)) {
            if (property(weightMap, (value(itE))) == property(label, sourceVertex(itE)) + property(label, targetVertex(itE))) {
                if (property(vertMap, sourceVertex(itE)) == true) addEdge(equalGraph, targetVertex(itE), sourceVertex(itE));
                else addEdge(equalGraph, sourceVertex(itE), targetVertex(itE));
            }
        }

        // Create a new matching
        clear(edges);
        matchSize = bipartiteMatching(equalGraph, vertMap, edges);
    }

    typedef typename Iterator<String<TEdges> >::Type TStringEdgeIter;
    TStringEdgeIter itSE = begin(edges);
    TStringEdgeIter itSEEnd = end(edges);
    TCargo sumWeight = 0;
    for(;itSE != itSEEnd; goNext(itSE)) sumWeight += property(weightMap, findEdge(g, (value(itSE)).i1, (value(itSE)).i2));

    return sumWeight;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TVertexMap, typename TWeightMap, typename TEdges>
typename Value<TWeightMap>::Type
weightedBipartiteMatching(String<TEdges> & edges,
                          Graph<TSpec> const & g,
                          TVertexMap const & vertMap,
                          TWeightMap const & weightMap)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Value<TWeightMap>::Type TCargo;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

    // Collect the two vertex sets, set1 is marked with false, set2 with true
    typedef String<TVertexDescriptor> TVertexSet;
    typedef typename Iterator<TVertexSet>::Type TVertexSetIter;
    TVertexSet set1;
    TVertexSet set2;
    TVertexIterator itV(g);
    for(;!atEnd(itV); goNext(itV)) {
        if (property(vertMap, value(itV)) == false) appendValue(set1, value(itV));
        else appendValue(set2, value(itV));
    }
    bool setIdentifier = true;      // Indicates what set needs more vertices
    TSize maxN = length(set1);
    if (maxN < length(set2)) {  maxN = length(set2); setIdentifier = false; }


    // Copy the original graph
    TGraph fullGraph;
    typedef typename EdgeType<TGraph>::Type TEdgeStump;
    resize(fullGraph.data_vertex, length(_getVertexString(g)), (TEdgeStump*) 0);
    fullGraph.data_id_managerV = g.data_id_managerV;
    TVertexMap myVertexMap = vertMap;
    resize(myVertexMap, maxN + maxN, setIdentifier);
    String<TCargo> myWeightMap;
    resize(myWeightMap, maxN * maxN, 0);
    TEdgeIterator itE(g);
    typedef std::pair<TVertexDescriptor, TVertexDescriptor> TEdge;
    typedef std::set<TEdge> TEdgeSet;
    TEdgeSet edgeSet;
    for(;!atEnd(itE); goNext(itE)) {
        TVertexDescriptor sV = sourceVertex(itE);
        TVertexDescriptor tV = targetVertex(itE);
        TEdgeDescriptor e = addEdge(fullGraph, sV, tV);
        if (sV < tV) edgeSet.insert(std::make_pair(sV, tV));
        else edgeSet.insert(std::make_pair(tV, sV));
        property(myWeightMap, e) = getProperty(weightMap, (value(itE)));
    }

    // Build a full graph
    if (setIdentifier == false) {
        TSize inc = maxN - length(set1);
        for(TSize i = 0; i< inc; ++i) appendValue(set1, addVertex(fullGraph));
    } else {
        TSize inc = maxN - length(set2);
        for(TSize i = 0; i<inc ; ++i) appendValue(set2, addVertex(fullGraph));
    }
    TVertexSetIter set1It = begin(set1);
    TVertexSetIter set1ItEnd = end(set1);
    for(;set1It != set1ItEnd; ++set1It) {
        TVertexSetIter set2It = begin(set2);
        TVertexSetIter set2ItEnd = end(set2);
        for(;set2It != set2ItEnd; ++set2It) {
            TVertexDescriptor sV = value(set1It);
            TVertexDescriptor tV = value(set2It);
            if (sV > tV) { TVertexDescriptor tmp = sV; sV = tV; tV = tmp; }
            if (edgeSet.find(std::make_pair(sV, tV)) == edgeSet.end()) addEdge(fullGraph, sV, tV);
        }
    }

    // Find a maximum weight matching
    String<TEdges> pseudo_edges;
    TCargo weight = _weightedBipartiteMatching(fullGraph, myVertexMap, myWeightMap, pseudo_edges);

    // Copy the relevant edges
    clear(edges);
    typedef typename Iterator<String<TEdges> >::Type TEdgeIter;
    TEdgeIter eIt = begin(pseudo_edges);
    TEdgeIter eItEnd = end(pseudo_edges);
    for(;eIt != eItEnd; ++eIt) {
        TVertexDescriptor sV = (value(eIt)).i1;
        TVertexDescriptor tV = (value(eIt)).i2;
        if (sV > tV) { TVertexDescriptor tmp = sV; sV = tV; tV = tmp; }
        if (edgeSet.find(std::make_pair(sV, tV)) != edgeSet.end()) appendValue(edges, value(eIt));
    }
    return weight;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_WEIGHTED_BIPARTITE_MATCHING_H_
