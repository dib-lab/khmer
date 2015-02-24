///A tutorial about finding Mums.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    // We begin with a StringSet that stores multiple strings.
    StringSet<String<char> > mySet;
    resize(mySet, 3);
    mySet[0] = "SeqAn is a library for sequence analysis.";
    mySet[1] = "The String class is the fundamental sequence type in SeqAn.";
    mySet[2] = "Subsequences can be handled with SeqAn's Segment class.";

    // Then we create an Index of this StringSet.
    typedef Index<StringSet<String<char> > > TMyIndex;
    TMyIndex myIndex(mySet);

    // To find maximal unique matches (Mums), we use the Mums Iterator and set the minimum MUM length to 3.
    Iterator<TMyIndex, Mums>::Type myMUMiterator(myIndex, 3);
    String<SAValue<TMyIndex>::Type> occs;

    while (!atEnd(myMUMiterator))
    {
        // A multiple match can be represented by the positions it occurs at in every sequence and its length.
        // getOccurrences@ returns an unordered sequence of pairs (seqNo,seqOfs) the match occurs at.
        occs = getOccurrences(myMUMiterator);
        //To order them ascending according seqNo we use orderOccurrences.
        orderOccurrences(occs);

        for (unsigned i = 0; i < length(occs); ++i)
            std::cout << getValueI2(occs[i]) << ", ";

        // repLength returns the length of the match.
        std::cout << repLength(myMUMiterator) << "   ";

        // The match string itself can be determined with representative.
        std::cout << "\t\"" << representative(myMUMiterator) << '\"' << std::endl;

        ++myMUMiterator;
    }

    return 0;
}
