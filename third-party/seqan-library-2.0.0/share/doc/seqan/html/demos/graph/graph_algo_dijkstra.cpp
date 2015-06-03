// A tutorial about the dijkstra's algorithm, once using an external map and once using an internal map.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;

    // Graph creation: 10 directed edges (0,1), (0,3), ...
    TGraph g;
    Size<TGraph>::Type numEdges = 10;
    VertexDescriptor<TGraph>::Type edges[] = {0, 1, 0, 3, 1, 2, 1, 3, 2, 4, 3, 1, 3, 2, 3, 4, 4, 0, 4, 2};
    addEdges(g, edges, numEdges);

    // One external property map: Weight map
    unsigned int weights[] =    {10, 5, 1, 2, 4, 3, 9, 2, 7, 6};
    String<unsigned int> weightMap;
    assignEdgeMap(weightMap, g, weights);

    // Out-parameters: Predecessor and distance map
    String<unsigned int> predMap;
    String<unsigned int> distMap;

    // Dijkstra from vertex 0 (single source shortest paths)
    dijkstra(predMap, distMap, g, 0, weightMap);

    // Output distances of shortest paths
    Iterator<TGraph, VertexIterator>::Type it(g);
    while (!atEnd(it))
    {
        std::cout << "Distance from 0 to " << getValue(it) << ": ";
        std::cout << getProperty(distMap, getValue(it)) << std::endl;
        goNext(it);
    }
    return 0;
}
