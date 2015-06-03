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
    Finder<TIndex> finder(index);

    std::cout << "The text has " << length(index) << " characters and consists of " << countSequences(index) <<
    " sequences." << std::endl;

    // The Finder object has a pointer to the first, current and last hit
    // Each consecutive call sets the current pointer to the appropriate hit
    while (find(finder, "MISS"))
        std::cout << "Hit at position: " << position(finder) << std::endl;

    return 0;
}
