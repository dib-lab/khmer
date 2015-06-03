#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef StringSet<String<char> >    TText;
    typedef Index<TText>                TIndex;
    typedef SAValue<TIndex>::Type       TSAValue;

    TText text;
    appendValue(text, "MISSISSIPPI");
    appendValue(text, "MYMISSISAHAPPY");

    TIndex index(text);
    Iterator<TIndex, TopDown<> >::Type it(index);

    goDown(it, "SSI");

    std::cout << "SSI occurs in " << getFrequency(it) << " sequences." << std::endl;

    // 1 alternative to print the location of hits
    String<TSAValue> saPositions = getOccurrences(it);
    for (unsigned i = 0; i < length(saPositions); ++i)
        std::cout << "Hit in sequence " << saPositions[i].i1 << " at position " << saPositions[i].i2 << std::endl;

    std::cout << "----------------------------" << std::endl;

    // 2 alternative to print the location of hits
    Pair<unsigned> hitInterval = range(it);
    for (; hitInterval.i1 < hitInterval.i2; ++hitInterval.i1)
        std::cout << "Hit in sequence " << getFibre(index, FibreSA())[hitInterval.i1].i1 <<
        " at position " << getFibre(index, FibreSA())[hitInterval.i1].i2 << std::endl;

    return 0;
}
