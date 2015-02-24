#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef Index<CharString> TIndex;

    TIndex index("MISSISSIPPI");
    Iterator<TIndex, TopDown<> >::Type it(index);

    goDown(it, "ISSI");
    std::cout << "The string " << representative(it) << " occurs " << range(it).i2 - range(it).i1 <<
    " times in MISSISSIPPI and has " << repLength(it) << " characters." << std::endl;

    // Note that goDown follows the path STARTING with a given text. It only stops at the next node. Therefore the
    // output for the following code is the same as above, even though the search string changed.
    goRoot(it);
    goDown(it, "ISS");

    std::cout << "The string " << representative(it) << " occurs " << range(it).i2 - range(it).i1 <<
    " times in MISSISSIPPI and has " << repLength(it) << " characters." << std::endl;

    return 0;
}
