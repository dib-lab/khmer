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

    // Fill external edge weight map.
    unsigned weights[] = {10, 5, 1, 2, 4, 3, 9, 2, 7, 6};
    String<unsigned> weightMap;
    assignEdgeMap(weightMap, g, weights);

    // Run Dijkstra's algorithm from vertex 0.
    String<unsigned> predMap;
    String<unsigned> distMap;
    dijkstra(predMap, distMap, g, 0, weightMap);

    // Print results to stdout.
    std::cout << "Single-Source Shortest Paths: \n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Path from 0 to " << getValue(it) << ": ";
        _printPath(g, predMap, (TVertexDescriptor) 0, getValue(it));
        std::cout << " (Distance: " << getProperty(distMap, getValue(it)) << ")\n";
        goNext(it);
    }

    // We can achieve the same thing using an internal map that is edge cargos.
    typedef unsigned int TEdgeCargo;
    typedef Directed<TEdgeCargo> TEdges;
    typedef Graph<TEdges> TCargoGraph;

    // Construct graph with the same edges as above.
    TCargoGraph cargoG;
    addEdges(cargoG, edges, numEdges);

    // Fill internal edge weight map.
    InternalPropertyMap<TEdgeCargo> intMap;
    assignEdgeMap(intMap, cargoG, weights);

    // Run Dijkstra's algorithm from vertex 0.
    clear(predMap);
    clear(distMap);
    dijkstra(predMap, distMap, cargoG, 0, intMap);

    // Print result to stdout.
    std::cout << "\nSingle-Source Shortest Paths: \n";
    typedef Iterator<TCargoGraph, VertexIterator>::Type TCargoVertexIterator;
    TCargoVertexIterator itC(cargoG);
    while (!atEnd(itC))
    {
        std::cout << "Path from 0 to " << getValue(itC) << ": ";
        _printPath(g, predMap, (TVertexDescriptor)0, getValue(itC));
        std::cout << " (Distance: " << getProperty(distMap, getValue(itC)) << ")\n";
        goNext(itC);
    }

    return 0;
}
