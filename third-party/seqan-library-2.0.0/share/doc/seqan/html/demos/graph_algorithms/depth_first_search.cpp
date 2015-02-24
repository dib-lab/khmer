#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create graph with 8 directed edges (0,3), (0,1), ...
    TSize numEdges = 8;
    TVertexDescriptor edges[] = {0, 3, 0, 1, 1, 4, 2, 4, 2, 5, 3, 1, 4, 3, 5, 5};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << "\n";

    // Create external property map for the vertex names and assign to graph.
    char names[] = {'u', 'v', 'w', 'x', 'y', 'z'};
    String<char> nameMap;
    assignVertexMap(nameMap, g, names);

    // Perform a DFS search.
    String<unsigned int> predMap;
    String<unsigned int> discoveryTimeMap;
    String<unsigned int> finishingTimeMap;
    depthFirstSearch(predMap, discoveryTimeMap, finishingTimeMap, g);

    // Write the result to stdout.
    std::cout << "Depth-First search: \n";
    typedef Iterator<Graph<>, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
        std::cout << "Discovery time = " << getProperty(discoveryTimeMap, getValue(it)) << ",";
        std::cout << "Finishing time = " << getProperty(finishingTimeMap, getValue(it)) << ",";
        typedef Value<String<unsigned int> >::Type TPredVal;
        TPredVal pre = getProperty(predMap, getValue(it));
        if (pre != getNil<TVertexDescriptor>())
            std::cout << "Predecessor = " << getProperty(nameMap, pre) << "\n";
        else
            std::cout << "Predecessor = nil" << "\n";
        goNext(it);
    }

    return 0;
}
