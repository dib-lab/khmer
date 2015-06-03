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

    // Fill external property map for each edge weight and vertex names.
    unsigned int weights[] = {4, 8, 8, 11, 7, 2, 4, 9, 14, 7, 6, 10, 1, 2};
    char names[] = {'a', 'b', 'c', 'd', 'i', 'e', 'h', 'g', 'f'};
    String<int> weightMap;
    assignEdgeMap(weightMap, g, weights);
    String<char> nameMap;
    assignVertexMap(nameMap, g, names);

    // Run Kruskal's algorithm.
    String<TVertexDescriptor> treeEdges;
    kruskalsAlgorithm(treeEdges, g, 0, weightMap);

    // Print the result to stdout.
    std::cout << "Minimum Spanning Tree (Kruskal's algorithm): \n"
              << "Tree Edges: ";
    typedef Iterator<String<TVertexDescriptor> >::Type TStrIterator;
    TStrIterator it = begin(treeEdges);
    TStrIterator itEnd = end(treeEdges);
    while (it != itEnd)
    {
        std::cout << "(" << getProperty(nameMap, getValue(it)) << ",";
        goNext(it);
        std::cout << getProperty(nameMap, getValue(it)) << "), ";
        goNext(it);
    }
    std::cout << "\n";

    return 0;
}
