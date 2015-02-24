#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Undirected<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create a graph with 10 undirected edges {0,1}, {0,4}, ...
    TSize numEdges = 10;
    TVertexDescriptor edges[] = {0, 1, 0, 4, 1, 5, 2, 5, 2, 6, 2, 3, 3, 6, 3, 7, 5, 6, 6, 7};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << "\n";

    // Create external property map for the vertex names and assign to graph.
    String<char> nameMap;
    char names[] = {'r', 's', 't', 'u', 'v', 'w', 'x', 'y'};
    assignVertexMap(nameMap, g, names);

    // Perform a BFS search starting from vertex with descriptor 1.
    String<unsigned int> predMap;
    String<unsigned int> distMap;
    breadthFirstSearch(predMap, distMap, g, 1);

    // Write the result to stdout.
    std::cout << "Breadth-First search: \n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
        if (getProperty(distMap, getValue(it)) == _getInfinityDistance(distMap))
            SEQAN_FAIL("Should never reach here!");
        else
            std::cout << "Level = " << getProperty(distMap, getValue(it));

        typedef Value<String<unsigned int> >::Type TPredVal;
        TPredVal pre = getProperty(predMap, getValue(it));
        if (pre != getNil<TVertexDescriptor>())
            std::cout << ", Predecessor = " << getProperty(nameMap, pre) << "\n";
        else
            std::cout << ", Predecessor = nil" << "\n";
        goNext(it);
    }

    return 0;
}
