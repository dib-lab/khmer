#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 10 directed edges (0,2), (0,1), ...
    TSize numEdges = 10;
    TVertexDescriptor edges[] = {0, 2, 0, 1, 1, 3, 1, 2, 2, 5, 2, 4, 2, 3, 3, 5, 3, 4, 4, 5};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph to stdout.
    std::cout << g << "\n";

    // Create external edge property map with edge weights.
    int weights[] = {3, 5, 6, 2, 2, 4, 7, 1, -1, -2};
    String<int> weightMap;
    assignEdgeMap(weightMap, g, weights);

    // Run DAG shortest path computation from vertex with descriptor 1.
    String<unsigned> predMap;
    String<unsigned> distMap;
    dagShortestPath(predMap, distMap, g, 1, weightMap);

    // Print result to stdout.
    std::cout << "Single-Source Shortest Paths in DAG: \n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Path from 1 to " << getValue(it) << ": ";
        _printPath(g, predMap, (TVertexDescriptor)1, getValue(it));
        std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")\n";
        goNext(it);
    }

    return 0;
}
