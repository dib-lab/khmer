#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef Index<CharString> TIndex;

    TIndex index("tobeornottobe");
    Finder<TIndex> finder(index);

    // The Finder object has a pointer to the first, current and last hit
    // Each consecutive call sets the current pointer to the appropriate hit
    while (find(finder, "to"))
        std::cout << "Hit at position: " << position(finder) << std::endl;

    return 0;
}
