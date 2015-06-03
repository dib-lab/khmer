#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 9 directed edges (0,3), (0,1), ...
    TSize numEdges = 9;
    TVertexDescriptor edges[] = {0, 3, 0, 1, 1, 2, 3, 2, 5, 7, 5, 6, 6, 7, 6, 3, 8, 7};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << "\n";

    // Create external property map with vertex labels.
    String<std::string> nameMap;
    std::string names[] = {"shirt", "tie", "jacket", "belt", "watch", "undershorts", "pants", "shoes", "socks"};
    assignVertexMap(nameMap, g, names);

    // Get vertex descriptor in topological sort order.
    String<TVertexDescriptor> order;
    topologicalSort(order, g);

    // Write the result to stdout.
    std::cout << "Topological sort: \n";
    typedef Iterator<String<TVertexDescriptor> >::Type TStringIterator;
    TStringIterator it = begin(order);
    TStringIterator itEnd = end(order);
    while (it != itEnd)
    {
        std::cout << getProperty(nameMap, getValue(it)) << ",";
        goNext(it);
    }
    std::cout << "\n";

    return 0;
}
