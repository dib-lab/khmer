#include <seqan/find.h>

using namespace seqan;

int main()
{
    CharString hstck = "I spy with my little eye something that is yellow";
    Finder<CharString> finder(hstck);

    Pattern<CharString, Horspool> p1("y");

    while (find(finder, p1))
        std::cout << "Hit at position: " << position(finder) << std::endl;

    goBegin(finder);    // move Finder to the beginning of the text
    clear(finder);      // reset Finder

    Pattern<CharString, Horspool> p2("t");
    while (find(finder, p2))
        std::cout << "Hit at position: " << position(finder) << std::endl;
}
