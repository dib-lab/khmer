//Examples for three seed extension algorithms.
#include <iostream>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
    String<char> a = "SEEDabXcdXefXXX";
    String<char> b = "SEEDabYcdefYYYY";

    Seed<Simple> seed1(0, 0, 4);          //left=0; length=4
    extendSeed(seed1, a, b, EXTEND_BOTH, MatchExtend());
    std::cout << "endPositionH(seed1) = " << endPositionH(seed1) << std::endl;
    std::cout << "endPositionV(seed1) = " << endPositionV(seed1) << std::endl;

    Seed<Simple> seed2(0, 0, 4);          //left=0; length=4
    Score<> scoring(1, -1, -1);
    extendSeed(seed2, a, b, EXTEND_BOTH, scoring, 2, UnGappedXDrop());
    std::cout << "endPositionH(seed2) = " << endPositionH(seed2) << std::endl;
    std::cout << "endPositionV(seed2) = " << endPositionV(seed2) << std::endl;

    Seed<Simple> seed3(0, 0, 4);          //left=0; length=4
    extendSeed(seed3, a, b, EXTEND_BOTH, scoring, 2, GappedXDrop());
    std::cout << "endPositionH(seed3) = " << endPositionH(seed3) << std::endl;
    std::cout << "endPositionV(seed3) = " << endPositionV(seed3) << std::endl;

    return 0;
}
