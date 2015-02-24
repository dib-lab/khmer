#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Undirected<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 14 undirected edges {0,1}, {0,6}, ...
    TSize numEdges = 14;
    TVertexDescriptor edges[] = {0, 1, 0, 6, 1, 2, 1, 6, 2, 3, 2, 4, 2, 8, 3, 5, 3, 8, 4, 6, 4, 7, 5, 8, 6, 7, 7, 8};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << std::endl;

    // Fill two external property maps for edge weights and vertex names.
    unsigned int weights[] = {4, 8, 8, 11, 7, 2, 4, 9, 14, 7, 6, 10, 1, 2};
    char names[] = {'a', 'b', 'c', 'd', 'i', 'e', 'h', 'g', 'f'};
    String<int> weightMap;
    assignEdgeMap(weightMap, g, weights);
    String<char> nameMap;
    assignVertexMap(nameMap, g, names);

    // Run Prim's algorithm.
    String<TVertexDescriptor> predMap;
    primsAlgorithm(predMap, g, 0, weightMap);

    // Print result to stdout.
    std::cout << "Minimum Spanning Tree (Prim's algorithm): \n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Path from " << getProperty(nameMap, 0) << " to "
                  << getProperty(nameMap, getValue(it)) << ": ";
        _printPath(g, predMap, (TVertexDescriptor)0, getValue(it), nameMap);
        std::cout << "\n";
        goNext(it);
    }

    return 0;
}
