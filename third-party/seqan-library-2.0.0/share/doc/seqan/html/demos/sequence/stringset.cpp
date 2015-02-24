#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    StringSet<String<char> > stringSet;
    appendValue(stringSet, "Hello World!");                    // Append string to the end of the string set.

    std::cout << "Number of elements: " << length(stringSet) << std::endl;

    resize(stringSet, 3, Exact());                             // Adapt the size of the container, while keeping existing values.

    String<char> seq = "To be or not to be!";
    stringSet[1] = seq;                                        // Use subscript operator to assign a new value.
    stringSet[2] = "A man, a plan, a canal - Panama!";

    std::cout << "Number of elements: " << length(stringSet) << std::endl;

    typedef Iterator<StringSet<String<char> >, Standard>::Type TIterator;
    for (TIterator it = begin(stringSet, Standard()); it != end(stringSet, Standard()); ++it)
        std::cout << "Element " << position(it, stringSet) << ": " << *it << std::endl;

    clear(stringSet);                                           // Clear the contents of the StringSet.

    std::cout << "Number of elements: " << length(stringSet) << std::endl;

    return 0;
}
