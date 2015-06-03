#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<char> myString = "tobeornottobe";

    typedef Index<String<char> > TIndex;
    TIndex index(myString);

//![iteration]
    Iterator<TIndex, TopDown<ParentLinks<Postorder> > >::Type it(index);

    goBegin(it);
    while (!atEnd(it))
    {
        std::cout << representative(it) << std::endl;
        ++it;
    }
//![iteration]

    return 0;
}
