#include <seqan/index.h>

using namespace seqan;

int main()
{
    typedef Index<CharString> TIndex;
    TIndex index("TATAA");

    Iterator<TIndex, TopDown<ParentLinks<> > >::Type itDefault;
    itDefault = begin(index, TopDown<ParentLinks<> >());

    while (!atEnd(itDefault))
    {
        std::cout << representative(itDefault) << std::endl;
        goNext(itDefault);
    }

    std::cout << "--------------------------------" << std::endl;

    Iterator<TIndex, TopDown<ParentLinks<Postorder> > >::Type itPostOrder;
    itPostOrder = begin(index, TopDown<ParentLinks<Postorder> >());

    while (!atEnd(itPostOrder))
    {
        std::cout << representative(itPostOrder) << std::endl;
        goNext(itPostOrder);
    }

    return 0;
}
