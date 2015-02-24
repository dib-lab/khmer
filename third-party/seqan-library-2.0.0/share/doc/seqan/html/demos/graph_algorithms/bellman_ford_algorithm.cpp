#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 10 directed edges (0,1), (0,3), ...
    TSize numEdges = 10;
    TVertexDescriptor edges[] = {0, 1, 0, 3, 1, 2, 1, 3, 2, 4, 3, 1, 3, 2, 3, 4, 4, 0, 4, 2};
    TGraph g;
    addEdges(g, edges, numEdges);

    // Print graph.
    std::cout << g << "\n";

    // Create external edge property map and assign to graph.
    unsigned weights[] =    {10, 5, 1, 2, 4, 3, 9, 2, 7, 6};
    String<unsigned> weightMap;
    assignEdgeMap(weightMap, g, weights);

    // Run Bellman-Ford algorithm from vertex 0.  NB: Ford-Fulkerson also
    // detects negative cycles.
    String<unsigned int> predMap;
    String<unsigned int> distMap;
    bool noNegativeCycle = bellmanFordAlgorithm(predMap, distMap, g, 0, weightMap);

    // Print result to stdout.
    std::cout << "Single-Source Shortest Paths: " << "\n"
              << "Graph without negative cycles? " << noNegativeCycle << "\n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Path from 0 to " << getValue(it) << ": ";
        _printPath(g, predMap, (TVertexDescriptor) 0, getValue(it));
        std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")\n";
        goNext(it);
    }

    return 0;
}
