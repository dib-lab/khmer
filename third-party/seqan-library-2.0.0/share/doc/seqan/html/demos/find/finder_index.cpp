#include <seqan/find.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    CharString hstck = "I spy with my little eye something that is yellow";
    Index<CharString, FMIndex<> > index(hstck);
    Finder<Index<CharString, FMIndex<> > > finder(hstck);

    while (find(finder, "y"))
        std::cout << "Hit at position: " << position(finder) << std::endl;

    clear(finder);      // reset Finder

    while (find(finder, "t"))
        std::cout << "Hit at position: " << position(finder) << std::endl;
}
