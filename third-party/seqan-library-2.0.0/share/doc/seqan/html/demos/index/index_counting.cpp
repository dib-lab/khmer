///An example to demonstrate  the functions countChildren and countOccurrences
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    //We begin with a String to store our sequence.
    String<char> myString = "How many wood would a woodchuck chuck. A woodchuck chucks as much wood as a woodchuck could";

    //Then we create an Index of this StringSet.
    typedef Index<String<char> > TMyIndex;
    TMyIndex myIndex(myString);

    // We will use a TopDown Iterator that supports parent links, ommits
    // empty edges and traverses the index in preorder to print out the number of
    // children at each node (not the number of leafs in the subtree).
    Iterator<TMyIndex, TopDown<ParentLinks<PreorderEmptyEdges> > >::Type tdIterator(myIndex);
    Size<TMyIndex>::Type count;

    while (!atEnd(tdIterator))
    {
        //We print out the representatives of all nodes that have more than 3
        //children and the number of occurrences. Also, we print a message if a node
        //is a leaf.
        count = countChildren(tdIterator);
        if (count >= 3)
        {
            std::cout << "Representative " << representative(tdIterator) << " has " <<  count << " children and ";
            std::cout << countOccurrences(tdIterator) << " occurrences " << std::endl;
        }
        if (isLeaf(tdIterator))
            std::cout << "The node is a leaf " << std::endl;

        tdIterator++;
    }

    return 0;
}
