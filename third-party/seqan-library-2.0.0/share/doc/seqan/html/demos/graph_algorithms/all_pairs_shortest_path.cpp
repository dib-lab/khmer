#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create a graph with 9 directed edges (0,1), (0,2), ...
    TSize numEdges = 9;
    TVertexDescriptor edges[] = {0, 1, 0, 2, 0, 4, 1, 3, 1, 4, 2, 1, 3, 0, 3, 2, 4, 3};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << std::endl;

    // Create a property map with edge weights.  Note that we can use negative
    // weights since the edges are directed and there are no cycles.
    int weights[] = {3, 8, -4, 1, 7, 4, 2, -5, 6};
    String<int> weightMap;
    assignEdgeMap(weightMap, g, weights);

    // Compute all-pairs shortest path.
    String<int> distMat;
    String<TVertexDescriptor> predMat;
    allPairsShortestPath(distMat, predMat, g, weightMap);

    // Print the result to stdout.
    unsigned int len = static_cast<unsigned>(std::sqrt((double)length(distMat)));
    for (TSize row = 0; row < len; ++row)
        for (TSize col = 0; col < len; ++col)
        {
            std::cout << row << "," << col << " (Distance="
                      << getValue(distMat, row * len + col) << "): ";
            _printAllPairsShortestPath(g, predMat, row, col);
            std::cout << "\n";
        }

    return 0;
}
