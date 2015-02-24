#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef Index<CharString> TIndex;

    TIndex index("tobeornottobe");
    Iterator<TIndex, TopDown<ParentLinks<> > >::Type it(index);

    do
    {
        // Print the letters from the root to the current node
        std::cout << representative(it) << std::endl;

        if (!goDown(it) && !goRight(it))
            while (goUp(it) && !goRight(it))
                ;
    }
    while (!isRoot(it));

    return 0;
}
