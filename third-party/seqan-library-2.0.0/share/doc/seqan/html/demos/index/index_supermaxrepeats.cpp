///A tutorial about finding supermaximal repeats.
//![includes]
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

//![includes]
//![init]
int main()
{
    String<char> myString = "How many wood would a woodchuck chuck.";

    typedef Index<String<char> > TMyIndex;
    TMyIndex myIndex(myString);

//![init]
//![iteration]
    Iterator<TMyIndex, SuperMaxRepeats>::Type myRepeatIterator(myIndex, 3);

    while (!atEnd(myRepeatIterator))
    {
        // A repeat can be represented by its length and positions it occurs at.
        // Function getOccurrences returns an unordered sequence of these positions
        // The length of this sequence, i.e. the repeat abundance can be obtained
        // from countOccurrences.
        for (unsigned i = 0; i < countOccurrences(myRepeatIterator); ++i)
            std::cout << getOccurrences(myRepeatIterator)[i] << ", ";

        // Function repLength returns the length of the repeat string.
        std::cout << repLength(myRepeatIterator) << "   ";

        // The repeat string itself can be determined with function representative.
        std::cout << "\t\"" << representative(myRepeatIterator) << '\"' << std::endl;

        ++myRepeatIterator;
    }

    return 0;
}

//![iteration]
