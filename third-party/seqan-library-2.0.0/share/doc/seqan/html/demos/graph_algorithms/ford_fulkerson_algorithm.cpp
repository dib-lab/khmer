#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 10 directed edges (0,1), (0,4), ...
    TSize numEdges = 10;
    TVertexDescriptor edges[] = {0, 1, 0, 4, 1, 2, 1, 4, 2, 3, 2, 4, 4, 1, 4, 5, 5, 2, 5, 3};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << "\n";

    // Create external property map for the edge capacities and assign to the graph.
    String<unsigned int> capMap;
    unsigned capacity[] =    {16, 13, 12, 10, 20, 9, 4, 14, 7, 4};
    assignEdgeMap(capMap, g, capacity);

    // Run the Ford-Fulkerson algorithm for maximum flow computation from source
    // vertex 0 to sink vertex 3.  valF is the value of the flow.
    String<unsigned int> flow;
    unsigned valF = fordFulkersonAlgorithm(flow, g, 0, 3, capMap);

    // Print the result to stdout.
    std::cout << "Ford-Fulkerson (Value of the flow = " << valF << ")\n";
    TEdgeIterator itEdge(g);
    for (; !atEnd(itEdge); goNext(itEdge))
        std::cout << "(" << sourceVertex(itEdge) << "," << targetVertex(itEdge) << "): "
                  << "Flow: " << getProperty(flow, getValue(itEdge)) << ", Capacity: "
                  << getProperty(capMap, getValue(itEdge)) << "\n";

    return 0;
}
