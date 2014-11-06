// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Graph - Algorithms
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Elementary graph algorithms
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Breadth-first search
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.breadthFirstSearch:
..cat:Graph
..summary:Implements a breadth-first search on a graph.
..remarks:Breadth-first search computes the distance from source to all reachable
vertices. It also produces a breath-first tree where each node has a predecessor / parent.
..signature:breadthFirstSearch(g, source, predecessor, distance)
..param.g:In-parameter:A graph.
...type:Spec.Undirected Graph
...type:Spec.Directed Graph
..param.source:In-parameter:A vertex descriptor.
...type:Metafunction.VertexDescriptor
...remarks:The breadth-first search is started from this vertex.
..param.predecessor:Out-parameter:A property map.
...remarks:The predecessor map stores implicitly the breadth-first tree.
..param.distance:Out-parameter:A property map.
...remarks:The distance map indicates at what depth a vertex was discovered.
..returns:void.
..see:Function.depthFirstSearch
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor, typename TPredecessorMap, typename TDistanceMap>
void
breadthFirstSearch(Graph<TSpec> const& g,
					 TVertexDescriptor const source,
					 TPredecessorMap& predecessor, 
					 TDistanceMap& distance)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TDistanceMap>::Type TDistVal;

	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);
	TPredVal nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
	TDistVal infDist = _getInfinityDistance(distance);
	
	String<bool> tokenMap;
	resizeVertexMap(g, tokenMap);
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(tokenMap, getValue(it), false);
		assignProperty(distance, getValue(it), infDist);
		assignProperty(predecessor, getValue(it), nilPred);
	}
	assignProperty(tokenMap, source, true);
	assignProperty(distance, source, 0);
	assignProperty(predecessor, source, nilPred);
	std::deque<TVertexDescriptor> queue;
	queue.push_back(source);
	
	// Bfs
	while (!queue.empty()) {
		TVertexDescriptor u = queue.front();
		queue.pop_front();
		typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
		TOutEdgeIterator itout(g,u);
		for(;!atEnd(itout);goNext(itout)) {
			TVertexDescriptor v = targetVertex(itout);
			if (getProperty(tokenMap, v) == false) {
				assignProperty(tokenMap, v, true);
				assignProperty(distance, v, getProperty(distance,u) + 1);
				assignProperty(predecessor, v, u);
				queue.push_back(v);
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Depth-first search
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor, typename TTokenMap, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap, typename TVal>
void
_dfsVisit(Graph<TSpec> const& g,
		   TVertexDescriptor const u,
		   TTokenMap& tokenMap,
		   TPredecessorMap& predecessor,
		   TDiscoveryTimeMap& disc,
		   TFinishingTimeMap& finish,
		   TVal& time)
{
	SEQAN_CHECKPOINT

	typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator;

	assignProperty(tokenMap, u, true);
	++time;
	assignProperty(disc, u, time);
	TAdjacencyIterator itad(g,u);
	for(;!atEnd(itad);goNext(itad)) {
		TVertexDescriptor v = getValue(itad);
		if (getProperty(tokenMap, v) == false) {
			assignProperty(predecessor, v, u);
			_dfsVisit(g, v, tokenMap, predecessor, disc, finish, time);
		}
	}
	++time;
	assignProperty(finish, u, time);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.depthFirstSearch:
..cat:Graph
..summary:Implements a depth-first search on a graph.
..remarks:In contrast to a breadth-first search the depth-first search is repeated from multiple sources if the graph is not connected.
Hence, depth-first search produces a depth-first forest. To ensure each vertex ends up in exactly one tree we need not just a distance but a
discovery and finishing time.
..signature:depthFirstSearch(g, predecessor, discovery, finish)
..param.g:In-parameter:A graph.
...type:Spec.Undirected Graph
...type:Spec.Directed Graph
..param.predecessor:Out-parameter:A property map.
...remarks:Predecessor subgraph produced by the depth-first search.
..param.discovery:Out-parameter:A property map.
...remarks:The discovery time of a vertex v.
..param.finish:Out-parameter:A property map.
...remarks:The time when v's adjacency list has been fully explored.
..returns:void.
..see:Function.breadthFirstSearch
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap>
void
depthFirstSearch(Graph<TSpec> const& g,
				   TPredecessorMap& predecessor,
				   TDiscoveryTimeMap& disc,
				   TFinishingTimeMap& finish)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPredecessorMap>::Type TPredVal;

	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,disc);
	resizeVertexMap(g,finish);
	TPredVal nilPred = getNil<TVertexDescriptor>();

	String<bool> tokenMap;
	resizeVertexMap(g, tokenMap);
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(tokenMap, getValue(it), false);
		assignProperty(predecessor, getValue(it), nilPred);
	}

	TSize time = 0;

	goBegin(it);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		if (getProperty(tokenMap, u) == false) {
			_dfsVisit(g, u, tokenMap, predecessor, disc, finish, time);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Topological sort
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.topologicalSort:
..cat:Graph
..summary:Performs a topological sort on a directed acyclic graph (DAG).
..remarks:A topological sort is a linear ordering of all its vertices such that if the graph contains an edge (u,v) then u appears before v in the ordering.
..signature:topologicalSort(g, topSort)
..param.g:In-parameter:A directed acyclic graph.
...type:Spec.Directed Graph
..param.topSort:Out-parameter:A topological ordering of the vertices.
...type:Class.String
..returns:void.
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor>
void
topologicalSort(Graph<TSpec> const & g,
				 String<TVertexDescriptor> & topSort)
{
	typedef typename Size<Graph<TSpec> >::Type TSize;

	// Variable definition.
	String<TSize> predMap;
	String<TSize> discoveryTimeMap;
	String<TSize> finishingTimeMap;
	
	// Perform DFS.
	depthFirstSearch(g, predMap, discoveryTimeMap, finishingTimeMap);
    SEQAN_ASSERT_EQ(numVertices(g), length(predMap));
    SEQAN_ASSERT_EQ(numVertices(g), length(discoveryTimeMap));
    SEQAN_ASSERT_EQ(numVertices(g), length(finishingTimeMap));

	// Order vertices.
	typedef ::std::pair<TSize, TVertexDescriptor> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	for (; !atEnd(it); goNext(it))
		q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));

	// Create topological order.
	resize(topSort, numVertices(g));
	TSize count = 0;
	while (!q.empty())
    {
		assignValue(topSort, count, q.top().second);
		q.pop();
		++count;
	}
    SEQAN_ASSERT_EQ(length(topSort), numVertices(g));
}

//////////////////////////////////////////////////////////////////////////////
// Strongly connected components
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.stronglyConnectedComponents:
..cat:Graph
..summary:Decomposes a directed graph into its strongly connected components.
..signature:stronglyConnectedComponents(g, components)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.components:Out-parameter:A property map.
...remarks:Each vertex is mapped to a component id. If two vertices share the same id they are in the same component.
..returns:$Size<TGraph>::Type$, number of strongly connected components.
..include:seqan/graph_algorithms.h
*/

template<typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
stronglyConnectedComponents(Graph<TSpec> const& g_source,
							  TComponents& components)
{
	// Initialization
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Value<TComponents>::Type TCompVal;
	resizeVertexMap(g_source,components);
	String<TSize> predMap;
	String<TSize> discoveryTimeMap;
	String<TSize> finishingTimeMap;

	// Dfs
	depthFirstSearch(g_source, predMap, discoveryTimeMap, finishingTimeMap);

	Graph<TSpec> g;
	transpose(g_source, g);

	// Second Dfs
	String<TSize> predecessor;
	String<TSize> disc;
	String<TSize> finish;
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,disc);
	resizeVertexMap(g,finish);
	TCompVal nilPred = getNil<TVertexDescriptor>();
	String<bool> tokenMap;
	resizeVertexMap(g, tokenMap);
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(components, getValue(it), nilPred);
		assignProperty(tokenMap, getValue(it), false);
		assignProperty(predecessor, getValue(it), nilPred);
	}

	// Order vertices
	typedef ::std::pair<TSize, TVertexDescriptor> TTimeVertexPair;
	std::priority_queue<TTimeVertexPair> q;
	goBegin(it);
	for(;!atEnd(it);++it) {
		q.push(std::make_pair(getProperty(finishingTimeMap, getValue(it)), getValue(it)));
	}

	TSize time = 0;
	TSize label = 0;
	while(!q.empty()) {
		TVertexDescriptor u = q.top().second;
		q.pop();
		if (getProperty(tokenMap, u) == false) {
			_dfsVisit(g, u, tokenMap, predecessor, disc, finish, time);
			TVertexIterator it_label(g);
			for(;!atEnd(it_label);goNext(it_label)) {
				if ((getProperty(tokenMap, getValue(it_label)) == true) &&
					(getProperty(components, getValue(it_label)) == nilPred)) {
					assignProperty(components, getValue(it_label), label);
				}
			}
			++label;
		}
	}

    return label;
}

//////////////////////////////////////////////////////////////////////////////
// Connected components
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor, typename TTokenMap, typename TComponents, typename TVal>
void
_connectedComponentVisit(Graph<TSpec> const& g,
		  TVertexDescriptor const u,
		  TTokenMap& tokenMap,
		  TComponents& components,
		  TVal& label)
{
	SEQAN_CHECKPOINT

	typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator;

	assignProperty(tokenMap, u, true);
	assignProperty(components, u, label);
	TAdjacencyIterator itad(g,u);
	for(;!atEnd(itad);goNext(itad)) {
		TVertexDescriptor v = getValue(itad);
		if (getProperty(tokenMap, v) == false) {
			_connectedComponentVisit(g, v, tokenMap, components, label);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.connectedComponents:
..cat:Graph
..summary:Decomposes an undirected graph into its connected components.
..signature:connectedComponents(g, components)
..param.g:In-parameter:An undirected graph.
...type:Spec.Undirected Graph
..param.components:Out-parameter:A property map.
...remarks:Each vertex is mapped to a component id. If two vertices share the same id they are in the same component.
..returns: The number of components.
..example:A simple example on how to use this function.
..example.code:
// Build Input.
Graph<Undirected<> > graph;
for (unsigned i = 0; i < 5; ++i)
    addVertex(graph);
addEdge(graph, 0, 1);
addEdge(graph, 0, 3);
addEdge(graph, 2, 4);
String<unsigned> components;
unsigned numComponents = 0;

// Call Algorithm.
numComponents = connectedComponents(g, components);

// Print Result.
std::cout << "Number of components: " << numComponents << std::endl;
std::cout << std::endl << "Vertex -> Component" << std::endl;
for (unsigned i = 0; i < length(components); ++i)
    std::cout << i << " -> " << components[i] << std::endl;
..example:The output now is:
..example.code:
Number of components: 2

Vertex -> Component
0 -> 0
1 -> 0
2 -> 1
3 -> 0
4 -> 1
..include:seqan/graph_algorithms.h
*/

template<typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
connectedComponents(Graph<TSpec> const& g_source,
					 TComponents& components)
{
	typedef typename Size<Graph<TSpec> >::Type TSize;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	clear(components);
	resizeVertexMap(g_source,components);
	
	// Initialization
	String<bool> tokenMap;
	resize(tokenMap, getIdUpperBound(_getVertexIdManager(g_source)), false);

	// Connected components
	TSize label = 0;
	TVertexIterator it(g_source);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		if (getProperty(tokenMap, u) == false) {
			_connectedComponentVisit(g_source, u, tokenMap, components, label);
			++label;
		}
	}
	return label;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Minimum Spanning Trees
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Prim's algorithm
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Function.primsAlgorithm:
..cat:Graph
..summary:Computes a minimum spanning tree on a graph.
..signature:primsAlgorithm(g, source, weight, predecessor)
..param.g:In-parameter:An undirected graph.
...type:Spec.Undirected Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:Edge weights.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a minimum spanning tree.
..returns:void.
..see:Function.kruskalsAlgorithm
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap>
void
primsAlgorithm(Graph<TSpec> const& g,
				TVertexDescriptor const source,
				TWeightMap const& weight,
				TPredecessorMap& predecessor)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TPredecessorMap>::Type TPred;
	typedef typename Value<TWeightMap>::Type TWeight;

	typedef ::std::pair<TWeight, TVertexDescriptor> TWeightVertexPair;
	std::priority_queue<TWeightVertexPair, std::vector<TWeightVertexPair>, std::greater<TWeightVertexPair> > q;
	
	// Initialization
	String<bool> tokenMap;
	String<TWeight> key;
	TPred nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
	TWeight infWeight = _getInfinityDistance(weight);
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,tokenMap);
	resizeVertexMap(g,key);

	TVertexIterator it(g);
	while(!atEnd(it)) {
		TVertexDescriptor u = getValue(it);
		if (u == source) q.push(std::make_pair(0, u));
		assignProperty(predecessor, u, nilPred);
		assignProperty(key, u, infWeight);
		assignProperty(tokenMap, u, false);
		goNext(it);
	}

	assignProperty(key, source, 0);
	while(!q.empty()) {
		TVertexDescriptor u = q.top().second;
		q.pop();
		if (getProperty(tokenMap, u)) continue;
		assignProperty(tokenMap, u, true);
		TOutEdgeIterator itOut(g,u);
		while(!atEnd(itOut)) {
			TVertexDescriptor v = targetVertex(itOut);
			TWeight w = getProperty(weight, getValue(itOut));
			if ((!getProperty(tokenMap, v)) &&
				(w < getProperty(key, v))) {
					assignProperty(predecessor, v, u);
					assignProperty(key, v, w);
					q.push(std::make_pair(w, v));
			}
			goNext(itOut);
		}
	}
}

template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap>
void
primsAlgorithmSpaceEfficient(Graph<TSpec> const& g,
							   TVertexDescriptor const source,
							   TWeightMap const& weight,
							   TPredecessorMap& predecessor)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TPredecessorMap>::Type TPred;
	typedef typename Value<TWeightMap>::Type TWeight;

	// Set-up the priority queue
	typedef Pair<TVertexDescriptor, TWeight> TKeyValue;
	typedef HeapTree<TKeyValue, std::less<TWeight>, KeyedHeap<> > TKeyedHeap;
	TKeyedHeap priorityQueue;
	
	// Initialization
	String<bool> tokenMap;
	TPred nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
	TWeight infWeight = _getInfinityDistance(weight);
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,tokenMap);

	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = value(it);
		heapInsert(priorityQueue, TKeyValue(u, infWeight));
		assignProperty(predecessor, u, nilPred);
		assignProperty(tokenMap, u, false);
	}
	heapChangeValue(priorityQueue, source, 0);

	// Iterate until queue is empty
	while(!empty(priorityQueue)) {
		TKeyValue kv = heapExtractRoot(priorityQueue);
		TVertexDescriptor u = kv.i1;
		assignProperty(tokenMap, u, true);
		if (kv.i2 == infWeight) continue;
		TOutEdgeIterator itOut(g,u);
		for(;!atEnd(itOut);goNext(itOut)) {
			TVertexDescriptor v = targetVertex(itOut);
			if (getProperty(tokenMap, v)) continue;
			TWeight w = getProperty(weight, getValue(itOut));
			if (w < heapGetValue(priorityQueue, v)) {
				assignProperty(predecessor, v, u);
				heapChangeValue(priorityQueue, v, w);
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Kruskal's algorithm
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template<typename TWeight, typename TPair>
struct LessPairI1_ :
	public ::std::unary_function<Pair<TWeight, TPair>, bool>
{
	inline bool 
	operator() (Pair<TWeight, TPair> const& a1, Pair<TWeight, TPair> const& a2) const {
		return (a1.i1 < a2.i1);
	}
};

/**
.Function.kruskalsAlgorithm:
..cat:Graph
..summary:Computes a minimum spanning tree on a graph.
..signature:kruskalsAlgorithm(g, source, weight, edges)
..param.g:In-parameter:An undirected graph.
...type:Spec.Undirected Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:Edge weights.
..param.edges:Out-parameter:Array of vertex descriptors.
...remarks:Array or string where two consecutive entries are an edge.
..returns:void.
..see:Function.primsAlgorithm
..include:seqan/graph_algorithms.h
*/

template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TEdges>
void
kruskalsAlgorithm(Graph<TSpec> const & g,
				   TVertexDescriptor const,
				   TWeightMap const & weight,
				   TEdges & edges)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TWeightMap>::Type TWeight;

	typedef Pair<TVertexDescriptor, TVertexDescriptor> TVertexPair;
	typedef Pair<TWeight, TVertexPair> TWeightEdgePair;
	typedef String<TWeightEdgePair>  TEdgeList;
	typedef typename Iterator<TEdgeList>::Type TEdgeListIter;
	TEdgeList edgeList;

	// Initialization
	reserve(edges, 2 * (numVertices(g) - 1));
    UnionFind<TVertexDescriptor> unionFind;
    resizeVertexMap(g, unionFind);
	
	// Sort the edges
	TEdgeIterator itE(g);
	for(;!atEnd(itE);goNext(itE)) appendValue(edgeList, TWeightEdgePair(getProperty(weight, getValue(itE)), TVertexPair(sourceVertex(itE),targetVertex(itE))));
	std::sort(begin(edgeList, Standard() ), end(edgeList, Standard() ), LessPairI1_<TWeight, TVertexPair>() );

	// Process each edge
	TEdgeListIter itEdgeList = begin(edgeList, Standard());
	TEdgeListIter itEdgeListEnd = end(edgeList, Standard());
	for (; itEdgeList != itEdgeListEnd; goNext(itEdgeList)) {
		TVertexDescriptor x = value(itEdgeList).i2.i1;
		TVertexDescriptor y = value(itEdgeList).i2.i2;

        if (findSet(unionFind, x) == findSet(unionFind, y))
            continue;

        appendValue(edges, x);
        appendValue(edges, y);
        joinSets(unionFind, findSet(unionFind, x), findSet(unionFind, y));
	}
}

// ----------------------------------------------------------------------------
// Weakly Connected Components
// ----------------------------------------------------------------------------

/**
.Function.weaklyConnectedComponents:
..cat:Graph
..summary:Compute weakly connected components of a directed graph.
..signature:weaklyConnectedComponents(g, components)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.components:Out-parameter:A property map.
...remarks:Each vertex is mapped to a component id. If two vertices share the same id they are in the same component.
..returns:$Size<TGraph>::Type$, number of weakly connected components.
..remarks:The running time is $O(n a(n, n))$ where $a$ is the inverse Ackermann function and thus almost linear. The union find data structure is used since the graph implementations do not allow the efficient iteration of in-edges.
..include:seqan/graph_algorithms.h
 */

template<typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
weaklyConnectedComponents(Graph<TSpec> const & g,
                          TComponents & components)
{
	SEQAN_CHECKPOINT;
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	// Initialization.
    UnionFind<TVertexDescriptor> unionFind;
    resizeVertexMap(g, unionFind);

    // Iterate over all edges, joining weakly connected components.
    for (TEdgeIterator itE(g); !atEnd(itE); goNext(itE))
        joinSets(unionFind, findSet(unionFind, sourceVertex(itE)), findSet(unionFind, targetVertex(itE)));

    // Count number of sets.
    TSize setCount = 0;
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV))
        setCount += (findSet(unionFind, *itV) == *itV);

    // Build a map from graph vertex descriptor to component id.
    TSize nextId = 0;
    clear(components);
    resizeVertexMap(g, components, setCount);  // setCount is sentinel value
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV)) {
        if (getProperty(components, findSet(unionFind, *itV)) == setCount)
            assignProperty(components, findSet(unionFind, *itV), nextId++);
        assignProperty(components, *itV, getProperty(components, findSet(unionFind, *itV)));
    }

    return setCount;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Single-Source Shortest Paths
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor, typename TNameMap>
inline void
_printPath(Graph<TSpec> const& g,
			TPredecessorMap const& predecessor,
			TVertexDescriptor const source,
			TVertexDescriptor const v,
			TNameMap const& nameMap)
{
	if (source == v) {
		std::cout << getProperty(nameMap, source);
	} else if (getProperty(predecessor, v) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
		std::cout << "No path from " << getProperty(nameMap, source) << " to " << getProperty(nameMap, v) << " exists.";
	} else {
		_printPath(g,predecessor, source, getProperty(predecessor, v), nameMap);
		std::cout << "," << getProperty(nameMap, v);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor1, typename TVertexDescriptor2>
inline void
_printPath(Graph<TSpec> const& g,
			TPredecessorMap const& predecessor,
			TVertexDescriptor1 const source,
			TVertexDescriptor2 const v)
{
	if (source == v) {
		std::cout << source;
	} else if (getProperty(predecessor, v) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
		std::cout << "No path from " << source << " to " << v << " exists.";
	} else {
		_printPath(g,predecessor, source, getProperty(predecessor, v));
		std::cout << "," << v;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor1, typename TVertexDescriptor2, typename TEdgeSet>
inline bool
_collectEdges(Graph<TSpec> const& g,
			   TPredecessorMap const& predecessor,
			   TVertexDescriptor1 const source,
			   TVertexDescriptor2 const v,
			   TEdgeSet& edgeSet)
{
	if ((TVertexDescriptor1) source == (TVertexDescriptor1) v) {
		return true;
	} else if (getProperty(predecessor, v) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
		return false;
	} else {
		edgeSet.insert(findEdge(g, getProperty(predecessor, v), v));
		return _collectEdges(g,predecessor, source, getProperty(predecessor, v), edgeSet);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor, typename TEdgeSet>
inline bool
_collectEdges(Graph<TSpec> const& g,
			   TPredecessorMap const& predecessor,
			   TVertexDescriptor const source,
			   TEdgeSet& edgeSet)
{
	typedef Iterator<Graph<Undirected<> >, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	for(;!atEnd(it); goNext(it)) {
		if (!_collectEdges(g, predecessor, source, value(it), edgeSet)) {
			edgeSet.clear();
			return false;
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
inline void 
_initializeSingleSource(Graph<TSpec> const& g,
						  TVertexDescriptor const source,
						  TWeightMap const& weight,
						  TPredecessorMap& predecessor, 
						  TDistanceMap& distance)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Value<TPredecessorMap>::Type TPredVal;
	typedef typename Value<TWeightMap>::Type TDistVal;
	TPredVal nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
	TDistVal infDist = _getInfinityDistance(weight);
	
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(distance, getValue(it), infDist);
		assignProperty(predecessor, getValue(it), nilPred);
	}
	assignProperty(distance, source, 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap, typename TVertexDescriptor, typename TEdgeDescriptor>
inline void 
_relax(Graph<TSpec> const& g,
	    TWeightMap const& weight,
		TPredecessorMap& predecessor, 
		TDistanceMap& distance,
		TVertexDescriptor const u,
		TEdgeDescriptor const e)
{
	TVertexDescriptor v = targetVertex(g,e);
	if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,e)) {
		assignProperty(distance, v, getProperty(distance,u) + getProperty(weight,e));
		assignProperty(predecessor, v, u);
	}
}



//////////////////////////////////////////////////////////////////////////////
// DAG Shortest Path
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.dagShortestPath:
..cat:Graph
..summary:Computes shortest paths from a single source in a directed acyclic graph (DAG).
..signature:dagShortestPath(g, source, weight, predecessor, distance)
..param.g:In-parameter:A directed acyclic graph.
...type:Spec.Directed Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:A weight map.
...remarks:In a directed acyclic graph edge weights can be negative because no cycles do exist.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a shortest-paths tree.
..param.distance:Out-parameter:A property map.
...remarks:Indicates for each vertex the distance from the source.
..returns:void.
..see:Function.bellmanFordAlgorithm
..see:Function.dijkstra
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void
dagShortestPath(Graph<TSpec> const& g,
				  TVertexDescriptor const source,
				  TWeightMap const& weight,
				  TPredecessorMap& predecessor,
				  TDistanceMap& distance)
{
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TStringIterator;

	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);

	// Topological sort
	String<TVertexDescriptor> order;
	topologicalSort(g, order);

	_initializeSingleSource(g, source, weight, predecessor, distance);

	//DAG Shortest Paths
	TStringIterator it = begin(order);
	while(!atEnd(it)) {
		TOutEdgeIterator itout(g, getValue(it));
		for(;!atEnd(itout);++itout) {
			_relax(g,weight,predecessor, distance, getValue(it), getValue(itout));
		}
		goNext(it);
	}
}


//////////////////////////////////////////////////////////////////////////////
// Bellman-Ford
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.bellmanFordAlgorithm:
..cat:Graph
..summary:Computes shortest paths from a single source in a directed graph.
..remarks:Edge weights may be negative in the Bellman-Ford algorithm.
The out parameters are only valid if the algorithm returns true.
..signature:bellmanFordAlgorithm(g, source, weight, predecessor, distance)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights may be negative.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a shortest-paths tree.
..param.distance:Out-parameter:A property map.
...remarks:Indicates for each vertex the distance from the source.
..returns:True if the graph has no negative weight cycles, false otherwise.
..see:Function.dagShortestPath
..see:Function.dijkstra
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
bool
bellmanFordAlgorithm(Graph<TSpec> const& g,
					   TVertexDescriptor const source,
					   TWeightMap const& weight,
					   TPredecessorMap& predecessor,
					   TDistanceMap& distance)
{
	SEQAN_CHECKPOINT
	typedef typename Size<Graph<TSpec> >::Type TSize;

	// Initialization
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);
	_initializeSingleSource(g, source, weight, predecessor, distance);

	// Run Bellman-Ford
	for(TSize i=0; i<numVertices(g) - 1; ++i) {
		TVertexIterator it(g);
		for(;!atEnd(it);goNext(it)) {
			TVertexDescriptor u = getValue(it);
			TOutEdgeIterator itout(g, u);
			for(;!atEnd(itout);++itout) {
				_relax(g,weight,predecessor, distance, u, getValue(itout));
			}
		}
	}

	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		TOutEdgeIterator itout(g, u);
		for(;!atEnd(itout);++itout) {
			TVertexDescriptor v = targetVertex(g, getValue(itout));
			if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,getValue(itout))) {
				return false;
			}	
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Dijkstra
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.dijkstra:
..cat:Graph
..summary:Computes shortest paths from a single source in a graph.
..remarks:Edge weights have to be nonnegative.
..signature:dijkstra(g, source, weight, predecessor, distance)
..param.g:In-parameter:A graph.
...type:Spec.Directed Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights have to be nonnegative.
..param.predecessor:Out-parameter:A property map.
...remarks:A property map that represents predecessor relationships among vertices. It determines a shortest-paths tree.
..param.distance:Out-parameter:A property map.
...remarks:Indicates for each vertex the distance from the source.
..returns:void
..see:Function.dagShortestPath
..see:Function.bellmanFordAlgorithm
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
void
dijkstra(Graph<TSpec> const& g,
		 TVertexDescriptor const source,
		 TWeightMap const& weight,
		 TPredecessorMap& predecessor,
		 TDistanceMap& distance)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Value<TDistanceMap>::Type TDistVal;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Initialization
	resizeVertexMap(g,predecessor);
	resizeVertexMap(g,distance);

	// S is initially empty
	String<bool> setS;
	resize(setS, getIdUpperBound(_getVertexIdManager(g)), false);

	// Set-up the priority queue
	typedef Pair<TVertexDescriptor, TDistVal> TKeyValue;
	typedef HeapTree<TKeyValue, std::less<TDistVal>, KeyedHeap<> > TKeyedHeap;
	TKeyedHeap priorityQueue;
	TDistVal infDist = _getInfinityDistance(weight);
	TVertexDescriptor nilVertex = getNil<typename VertexDescriptor<TGraph>::Type>();
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(predecessor, value(it), nilVertex);
		assignProperty(distance, value(it), infDist);
		heapInsert(priorityQueue, TKeyValue(value(it), infDist));
	}
	assignProperty(distance, source, 0);
	heapChangeValue(priorityQueue, source, 0);

	// Run Dijkstra
	while (!empty(priorityQueue)) {
		// Extract min
		TVertexDescriptor u = heapExtractRoot(priorityQueue).i1;
		assignProperty(setS, u, true);
		TOutEdgeIterator itout(g, u);
		for(;!atEnd(itout);++itout) {
			TVertexDescriptor v = targetVertex(itout);
			if (property(setS, v) == true) continue;
			if (getProperty(distance, v) > getProperty(distance,u) + getProperty(weight,value(itout))) {
				assignProperty(distance, v, getProperty(distance,u) + getProperty(weight,value(itout)));
				assignProperty(predecessor, v, u);
				heapChangeValue(priorityQueue, v, getProperty(distance,u) + getProperty(weight,value(itout)));
			}
		}
	}
}

//template<typename TSpec, typename TVertexDescriptor, typename TWeightMap, typename TPredecessorMap, typename TDistanceMap>
//void 
//dijkstra(Graph<TSpec> const& g,
//		 TVertexDescriptor const source,
//		 TWeightMap const& weight,
//		 TPredecessorMap& predecessor, 
//		 TDistanceMap& distance)
//{
//	SEQAN_CHECKPOINT
//	typedef Graph<TSpec> TGraph;
//	typedef typename Size<TGraph>::Type TSize;
//	typedef typename Value<TDistanceMap>::Type TDistVal;
//	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
//	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
//	
//	// Initialization
//	resizeVertexMap(g,predecessor);
//	resizeVertexMap(g,distance);
//
//	_initializeSingleSource(g, source, weight, predecessor, distance);
//	
//	String<bool> setS;
//	resizeVertexMap(g, setS);
//	TVertexIterator it(g);
//	for(;!atEnd(it);++it) {
//		assignProperty(setS, getValue(it), false);
//	}
//	TDistVal infDist = _getInfinityDistance(weight);
//	TVertexDescriptor nilVertex = getNil<typename VertexDescriptor<TGraph>::Type>();
//
//	// Run Dijkstra
//	TSize count = numVertices(g);
//	while (count > 0) {
//		// Extract min
//		TDistVal min = infDist;
//		TVertexDescriptor u = nilVertex;
//		TVertexIterator it_find(g);
//		for(;!atEnd(it_find);++it_find) {
//			if(getProperty(setS,getValue(it_find))==true) continue;
//			if ((u == nilVertex) ||
//				(getProperty(distance,getValue(it_find))<getProperty(distance,u))) {
//					u = getValue(it_find);
//					min = getProperty(distance,getValue(it_find));
//			}
//		}
//		assignProperty(setS, u, true);
//		TOutEdgeIterator itout(g, u);
//		for(;!atEnd(itout);++itout) {
//			_relax(g,weight,predecessor, distance, u, getValue(itout));
//		}
//		--count;
//	}
//}





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// All-Pairs shortest paths
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessor, typename TVertexDescriptor>
inline void
_printAllPairsShortestPath(Graph<TSpec> const& g,
							   TPredecessor& predecessor, 
							   TVertexDescriptor const i,
							   TVertexDescriptor const j)
{
	typedef typename Size<TPredecessor>::Type TSize;
	TSize len = getIdUpperBound(g.data_id_managerV);
	if (i==j) {
		std::cout << i;
	} else if (getValue(predecessor, i*len+j) == getNil<typename VertexDescriptor<Graph<TSpec> >::Type>()) {
		std::cout << "No path from " << i << " to " << j << " exists.";
	} else {
		_printAllPairsShortestPath(g,predecessor, i, (TVertexDescriptor) getValue(predecessor, i*len+j));
		std::cout << "," << j;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
_initializeAllPairs(Graph<TSpec> const& g,
						TWeightMap const& weight,
						TMatrix& matrix,
						TPredecessor& predecessor)
{
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TWeightMap>::Type TWeightVal;
	typedef typename Value<TPredecessor>::Type TPredVal;
	
	// Create adjacency-like matrix
	TSize len = getIdUpperBound(g.data_id_managerV);
	resize(matrix, len * len);
	resize(predecessor, len * len);
	TWeightVal infWeight = _getInfinityDistance(weight);
	TPredVal nilPred = getNil<TVertexDescriptor>();
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			if (row != col) assignValue(matrix, row*len + col, infWeight);
			else assignValue(matrix, row*len + col, 0);
			assignValue(predecessor, row*len + col, nilPred);
		}
	}

	// Include edge weights and initial predecessors
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor u = getValue(it);
		TOutEdgeIterator itout(g, u);
		for(;!atEnd(itout);++itout) {
			TVertexDescriptor v = targetVertex(g,getValue(itout));
			assignValue(matrix, u*len + v, getProperty(weight, getValue(itout)));
			assignValue(predecessor, u*len + v, u);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TPredecessor, typename TInfDist>
void 
_extendShortestPaths(TMatrix& local,
					   TMatrix& w,
					   TPredecessor& predecessor,
					   TInfDist const infDist)
{
	typedef typename Value<TMatrix>::Type TMatrixVal;
	typedef typename Value<TPredecessor>::Type TPredVal;
	typedef typename Size<TMatrix>::Type TSize;
	TMatrix oldLocal = local;
	TPredecessor oldPredecessor = predecessor;
	TSize len = (TSize) std::sqrt((double) length(oldLocal));
	for(TSize i = 0; i<len;++i) {
		for(TSize j = 0; j<len;++j) {
			if (i==j) continue;
			assignValue(local, i*len+j,infDist);
			TPredVal ind = 0;
			for(TSize k = 0; k<len;++k) {
				TMatrixVal min1 = getValue(local, i*len+j);
				TMatrixVal min2 = getValue(oldLocal, i*len+k) + getValue(w, k*len + j);
				if (min2 < min1) {
					assignValue(local, i*len+j,min2);
					ind = k;
				}
			}
			if (getValue(oldLocal, i*len+j) > getValue(local, i*len+j)) {
				assignValue(predecessor, i*len+j,ind);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// All-Pairs shortest path
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Function.allPairsShortestPath:
..cat:Graph
..summary:Finds shortest paths between all pairs of vertices in a graph.
..signature:allPairsShortestPath(g, weight, distance, predecessor)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights may be negative.
..param.distance:Out-parameter:A matrix with distances.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the distance from vertex i to vertex j.
..param.predecessor:Out-parameter:A matrix with predecessors.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the predecessor of j on a shortest path from vertex i to vertex j.
You can use _printAllPairsShortestPath(g, predecessor, i, j) to print the shortest path from i to j.
..returns:void
..see:Function.floydWarshallAlgorithm
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
allPairsShortestPath(Graph<TSpec> const& g,
						TWeightMap const& weight,
						TMatrix& distMatrix,
						TPredecessor& predecessor)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TWeightMap>::Type TWeightVal;
	TWeightVal infWeight = _getInfinityDistance(weight);

	// Initialize first distance matrix
	_initializeAllPairs(g,weight,distMatrix,predecessor);

	TSize len = (TSize) sqrt((double) length(distMatrix));
	TMatrix local = distMatrix;
	for(TSize m=2;m<len;++m) {
		_extendShortestPaths(local,distMatrix,predecessor, infWeight);
	}
	distMatrix = local;
}


//////////////////////////////////////////////////////////////////////////////
// Floyd-Warshall
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Function.floydWarshallAlgorithm:
..cat:Graph
..summary:Finds shortest paths between all pairs of vertices in a graph.
..signature:floydWarshallAlgorithm(g, weight, distance, predecessor)
..remarks:The graph must be free of negative-weight cycles.
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.weight:In-parameter:A weight map.
...remarks:A property map with edge weights. Edge weights may be negative.
..param.distance:Out-parameter:A matrix with distances.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the distance from vertex i to vertex j.
..param.predecessor:Out-parameter:A matrix with predecessors.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates the predecessor of j on a shortest path from vertex i to vertex j.
You can use _printAllPairsShortestPath(g, predecessor, i, j) to print the shortest path from i to j.
..returns:void
..see:Function.allPairsShortestPath
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TWeightMap, typename TMatrix, typename TPredecessor>
void 
floydWarshallAlgorithm(Graph<TSpec> const& g,
			   TWeightMap const& weight,
			   TMatrix& distMatrix,
			   TPredecessor& predecessor)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TMatrix>::Type TSize;
	typedef typename Value<TMatrix>::Type TMatrixVal;

	// Initialize first distance matrix
	_initializeAllPairs(g,weight,distMatrix,predecessor);

	// Floyd-Warshall
	TSize len = (TSize) std::sqrt((double) length(distMatrix));
	TMatrix local = distMatrix;
	for(TSize k=0;k<len;++k) {
		for(TSize i=0;i<len;++i) {
			for(TSize j=0;j<len;++j) {
				TMatrixVal min1 = getValue(distMatrix, i*len+j);
				TMatrixVal min2 = getValue(distMatrix, i*len+k) + getValue(distMatrix, k*len + j);
				if (min2 < min1) {
					assignValue(local, i*len+j,min2);
					assignValue(predecessor, i*len+j,getValue(predecessor, k*len+j));
				} else {
					assignValue(local, i*len+j,min1);
					assignValue(predecessor, i*len+j, getValue(predecessor, i*len+j));
				}
			}
		}
		distMatrix=local;
	}
}

//////////////////////////////////////////////////////////////////////////////
// Transitive Closure
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/**
.Function.transitiveClosure:
..cat:Graph
..summary:Determines whether there is a path between any two given vertices or not.
..signature:transitiveClosure(g, closure)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.closure:Out-parameter:A matrix which indicates the closure.
...type:Class.Matrix
...remarks:Entry (i,j) in this matrix indicates whether there is a path from i to j in the graph or not.
..returns:void
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TMatrix>
void
transitiveClosure(Graph<TSpec> const& g,
				   TMatrix& closure)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize first closure matrix
	getAdjacencyMatrix(g,closure);
	TSize len = (TSize) std::sqrt((double) length(closure));
	for (TSize diag=0;diag < len;++diag) assignValue(closure, diag*len+diag,1);

	// Transitive Closure
	TMatrix local = closure;
	for (TSize k=0;k<len;++k) {
		for(TSize i=0;i<len;++i) {
			for(TSize j=0;j<len;++j) {
				bool t_ij = static_cast<int>(getValue(closure, i*len+j)) > 0;
                bool t_ik = static_cast<int>(getValue(closure, i*len+k)) > 0;
                bool t_kj = static_cast<int>(getValue(closure, k*len+j)) > 0;
				assignValue(local, i*len+j, t_ij || (t_ik && t_kj));
			}
		}
		closure = local;
	}
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Maximum Flow
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
template<typename TSpec, typename TCapMap, typename TFlowMap, typename TResidualGraph>
void
_buildResidualGraph(Graph<TSpec> const& g,
					  TCapMap const& capacity,
					  TFlowMap const& flow,
					  TResidualGraph& rG)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TFlowMap>::Type TFlow;
	typedef typename Value<TCapMap>::Type TCap;

	clear(rG);
	TVertexIterator itV(g);
	for(;!atEnd(itV);goNext(itV)) {
		_createVertices(rG, getValue(itV));
	}

	TEdgeIterator itE(g);
	for(;!atEnd(itE);goNext(itE)) {
		typedef typename EdgeDescriptor<TResidualGraph>::Type TEdgeDescriptor;
		TFlow f = getProperty(flow, getValue(itE));
		TCap cap = getProperty(capacity, getValue(itE));
		if (f > 0) {
			TEdgeDescriptor e_rG = findEdge(rG, targetVertex(itE), sourceVertex(itE));
			if (e_rG == 0) addEdge(rG, targetVertex(itE), sourceVertex(itE), f);
			else cargo(e_rG) += f;
		}
		if (f < cap) {
			TEdgeDescriptor e_rG = findEdge(rG, sourceVertex(itE), targetVertex(itE));
			if (e_rG == 0) addEdge(rG, sourceVertex(itE), targetVertex(itE), cap - f);
			else cargo(e_rG) += cap - f;			
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPredecessorMap, typename TVertexDescriptor>
inline typename Size<Graph<TSpec> >::Type
_getMinimumAug(Graph<TSpec> const& rG,
				 TPredecessorMap& predecessor,
				 TVertexDescriptor const source,
				 TVertexDescriptor sink)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef TSize TFlow;
	typedef typename Iterator<String<TVertexDescriptor>, Rooted>::Type TIterator;
	
	// Build secondary predecessor map just containing the path
	TVertexDescriptor nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
	String<TVertexDescriptor> predMap;
	resizeVertexMap(rG, predMap);
	TIterator it = begin(predMap);
	for(;!atEnd(it);goNext(it)) {
		*it = nilPred;
	}

	// Find minimum flow
	TVertexDescriptor pred = getProperty(predecessor, sink);
	TFlow f = getCargo(findEdge(rG, pred,sink));
	assignProperty(predMap, sink, pred);
	while(pred != source) {
		sink = pred;
		pred = getProperty(predecessor, sink);
		TFlow f2 = getCargo(findEdge(rG, pred,sink));
		assignProperty(predMap, sink, pred);
		if (f2 < f) f = f2;
	}

	// Just return the augmenting path
	predecessor = predMap;
	return f;
}

//////////////////////////////////////////////////////////////////////////////
// Ford Fulkerson
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Function.fordFulkersonAlgorithm:
..cat:Graph
..summary:Computes a maximum flow in a directed graph.
..signature:fordFulkersonAlgorithm(g, source, sink, capacity, flow)
..param.g:In-parameter:A directed graph.
...type:Spec.Directed Graph
..param.source:In-parameter:A source vertex.
...type:Metafunction.VertexDescriptor
..param.sink:In-parameter:A sink vertex.
...type:Metafunction.VertexDescriptor
..param.capacity:In-parameter:A property map of edge capacities.
..param.flow:Out-parameter:A property map with the flow of each edge.
..returns:The value of the flow.
..include:seqan/graph_algorithms.h
*/
template<typename TSpec, typename TVertexDescriptor, typename TCapMap, typename TFlowMap>
typename Value<TFlowMap>::Type
fordFulkersonAlgorithm(Graph<TSpec> const& g,
			   TVertexDescriptor const source,
			   TVertexDescriptor const sink,
			   TCapMap const& capacity,
			   TFlowMap& flow)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TFlowMap>::Type TFlow;

	// Initialization
	TVertexDescriptor nilPred = getNil<typename VertexDescriptor<TGraph>::Type>();
	resizeEdgeMap(g,flow);
	TEdgeIterator itE(g);
	for(;!atEnd(itE);goNext(itE)) {
		assignProperty(flow, getValue(itE), 0);
	}

	// Build the residual graph
	Graph<Directed<TFlow> > rG;
	_buildResidualGraph(g,capacity, flow, rG);

		
	// Determine whether the sink is reachable
	String<TVertexDescriptor> predMap;
	String<TVertexDescriptor> distMap;
	breadthFirstSearch(rG, source, predMap, distMap);
	
	while (getProperty(predMap, sink) != nilPred) {
		TFlow inc = _getMinimumAug(rG, predMap, source, sink);
		TEdgeIterator itEdge(g);
		for(;!atEnd(itEdge);goNext(itEdge)) {
			TVertexDescriptor u = sourceVertex(itEdge);
			TVertexDescriptor v = targetVertex(itEdge);
			TEdgeDescriptor e = getValue(itEdge);
			if (getProperty(predMap, v) == u) assignProperty(flow, e, getProperty(flow, e) + inc);
			if (getProperty(predMap, u) == v) assignProperty(flow, e, getProperty(flow, e) - inc);
		}
		// Build the residual graph
		_buildResidualGraph(g,capacity, flow, rG);
		// Determine whether the sink is reachable
		clear(predMap);
		clear(distMap);
		breadthFirstSearch(rG, source, predMap, distMap);
	}

	TFlow valF = 0;
	TOutEdgeIterator itOutEdge(g, source);
	for(;!atEnd(itOutEdge);goNext(itOutEdge)) {
		valF += getProperty(flow, getValue(itOutEdge));
	}
	return valF;
}

















//////////////////////////////////////////////////////////////////////////////
// ToDo: Not yet tested, use with care
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// Matching
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Path Growing Algorithm
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TWeightMap, typename TEdgeMap>
typename Value<TWeightMap>::Type
pathGrowingAlgorithm(Graph<TSpec>& g,
					   TWeightMap const& weightMap,
					   TEdgeMap& edgeMap1)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename Value<TWeightMap>::Type TValue;
	typedef typename Size<Graph<TSpec> >::Type TSize;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Make a copy of the graph
	TGraph mutant(g);

	// Initialy not a single edge is selected
	resize(edgeMap1, getIdUpperBound(_getEdgeIdManager(g)), false);
	TEdgeMap edgeMap2 = edgeMap1;
	TValue edgeMap1Sum = 0;
	TValue edgeMap2Sum = 0;
	
	// Run the algorithm
	TSize i = 1;
	while (numEdges(mutant) > 0) {
		TVertexIterator itVert(mutant);
		while (outDegree(mutant, *itVert) < 1) goNext(itVert);
		TVertexDescriptor x = *itVert;
		TVertexDescriptor y;
		while (outDegree(mutant, x) >= 1) {
			TOutEdgeIterator itOut(mutant, x);
			TEdgeDescriptor e = *itOut;
			TValue max = getProperty(weightMap, e);
			y = targetVertex(itOut);
			goNext(itOut);
			for(;!atEnd(itOut);++itOut) {
				if (getProperty(weightMap, *itOut) > max) {
					e = *itOut;
					max = getProperty(weightMap, e);
					y = targetVertex(itOut);
				}
			}
			if (i == 1) {
				// Mark the edge for m1
				assignProperty(edgeMap1, e, true);
				edgeMap1Sum += max;
			} else {
				// Mark the edge for m2
				assignProperty(edgeMap2, e, true);
				edgeMap2Sum += max;
			}
			i = 3 - i;
			removeVertex(mutant, x);
			x = y;
		}
	}
	

	// Check whether we have to swap bool arrays
	if (edgeMap2Sum > edgeMap1Sum) {
		edgeMap1Sum = edgeMap2Sum;
		edgeMap1 = edgeMap2;
	}

	return edgeMap1Sum;
}




//////////////////////////////////////////////////////////////////////////////
// Weighted bipartite Matching
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexMap, typename TEdges>
inline typename Size<Graph<TSpec> >::Type
bipartiteMatching(Graph<TSpec>& g,			  
					TVertexMap& vertMap,
					String<TEdges>& edges)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIter;

	clear(edges);
	TVertexDescriptor source = addVertex(g);
	TVertexDescriptor target = addVertex(g);
	TVertexIter itV(g);
	for(;!atEnd(itV); goNext(itV)) {
		if ((value(itV) != source) && (value(itV) != target)) {
			if (getProperty(vertMap, value(itV)) == false) {
				addEdge(g, source, value(itV));
			} else {
				addEdge(g, value(itV), target);
			}
		}
	}

	// Use Ford-Fulkerson to determine a matching
	String<TSize> capMap;	
	resizeEdgeMap(g,capMap);
	typedef typename Iterator<String<TSize> >::Type TCapIter;
	TCapIter capIt = begin(capMap);
	TCapIter capItEnd = end(capMap);
	for(;capIt != capItEnd; ++capIt) value(capIt) = 1;
	String<TSize> flow;	
	TSize valF = fordFulkersonAlgorithm(g, source, target, capMap, flow);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(g);
	for(;!atEnd(itEdge);goNext(itEdge)) {
		if (getProperty(flow, getValue(itEdge)) == 1) {
			TVertexDescriptor sV = sourceVertex(itEdge);
			TVertexDescriptor tV = targetVertex(itEdge);
			if ((sV != source) && (tV != target)) appendValue(edges, TEdges(sV, tV));
		}
	}
	removeVertex(g, source);
	removeVertex(g, target);

	return valF;
}


/////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexMap, typename TWeightMap, typename TEdges>
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
	resizeVertexMap(g, label);
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
		resizeVertexMap(g, reverseMatchMap);
		
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
			if (property(vertMap, sV) == true) {	TVertexDescriptor tmp = sV;	sV = tV; tV = tmp;	}
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

template<typename TSpec, typename TVertexMap, typename TWeightMap, typename TEdges>
inline typename Value<TWeightMap>::Type
weightedBipartiteMatching(Graph<TSpec>& g,
							TVertexMap& vertMap,
							TWeightMap& weightMap,
							String<TEdges>& edges)
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
	bool setIdentifier = true;		// Indicates what set needs more vertices
	TSize maxN = length(set1);
	if (maxN < length(set2)) {	maxN = length(set2); setIdentifier = false; }


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


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
