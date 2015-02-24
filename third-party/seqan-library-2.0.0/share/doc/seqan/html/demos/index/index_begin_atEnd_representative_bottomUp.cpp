#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef Index<CharString> TIndex;
    TIndex index("TATAA");

    Iterator<TIndex, BottomUp<Postorder> >::Type itDefault;
    itDefault = begin(index, BottomUp<Postorder>());

    while (!isRoot(itDefault))
    {
        std::cout << representative(itDefault) << std::endl;
        goNext(itDefault);
    }

    return 0;
}
