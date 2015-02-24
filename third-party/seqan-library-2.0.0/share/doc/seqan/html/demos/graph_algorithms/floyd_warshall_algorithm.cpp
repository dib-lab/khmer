#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 9 directed edges (0,1), (0,2)
    TSize numEdges = 9;
    TVertexDescriptor edges[] = {0, 1, 0, 2, 0, 4, 1, 3, 1, 4, 2, 1, 3, 0, 3, 2, 4, 3};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << "\n";

    // Fill external property map with edge weights and assign to graph.
    int weights[] = {3, 8, -4, 1, 7, 4, 2, -5, 6};
    String<int> weightMap;
    assignEdgeMap(weightMap, g, weights);

    // Run Floyd-Warshall algorithm.
    String<int> distMat;
    String<TVertexDescriptor> predMat;
    floydWarshallAlgorithm(distMat, predMat, g, weightMap);

    // Print result to stdout.
    unsigned int len = static_cast<unsigned>(std::sqrt((double)length(distMat)));
    for (TSize row = 0; row < len; ++row)
        for (TSize col = 0; col < len; ++col)
        {
            std::cout << row << "," << col << " (Distance="
                      << getValue(distMat, row * len + col) << "): ";
            _printAllPairsShortestPath(g, predMat, row, col);
            std::cout << std::endl;
        }

    return 0;
}
