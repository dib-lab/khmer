#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Size<TGraph>::Type TSize;

    // Create grap with 14 directed edges (1,0), (0,4), ...
    TSize numEdges = 14;
    TVertexDescriptor edges[] = {1, 0, 0, 4, 2, 1, 4, 1, 5, 1, 6, 2, 3, 2, 2, 3, 7, 3, 5, 4, 6, 5, 5, 6, 7, 6, 7, 7};
    TGraph g;
    addEdges(g, edges, numEdges);
    // Print graph.
    std::cout << g << std::endl;

    // Create external property map with vertex names and assign to graph.
    String<char> nameMap;
    char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    assignVertexMap(nameMap, g, names);

    // Compute strongly connected components.
    String<unsigned int> component;
    stronglyConnectedComponents(component, g);

    // Print result to stdout.
    std::cout << "Strongly Connected Components: \n";
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    while (!atEnd(it))
    {
        std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": \n"
                  << "Component = " << getProperty(component, getValue(it)) << "\n";
        goNext(it);
    }

    return 0;
}
