    #include <iostream>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed> TSeedSet;
    typedef Iterator<TSeedSet>::Type TIterator;

    DnaString seqH = "ggatACGTccTTCGAtACccTGGTg";
    DnaString seqV = "ttccgACGTgTTCGAgACtgaGGTca";

    TSeed seed0(4, 5, 4);
    TSeed seed1(10, 10, 5);
    TSeed seed2(14, 14, 4);
    TSeed seed3(21, 21, 3);

    TSeedSet seedSet;

    addSeed(seedSet, seed0, Single());
    addSeed(seedSet, seed1, Single());
    addSeed(seedSet, seed2, Single());
    addSeed(seedSet, seed3, Single());

    std::cout << "Single Method:" << std::endl;
    for (TIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
        std::cout << "Seed: " << *it << std::endl;
    std::cout << std::endl;

    clear(seedSet);
    if (!addSeed(seedSet, seed0, 2, Merge()))
        addSeed(seedSet, seed0, Single());
    if (!addSeed(seedSet, seed1, 2, Merge()))
        addSeed(seedSet, seed1, Single());
    if (!addSeed(seedSet, seed2, 2, Merge()))
        addSeed(seedSet, seed2, Single());
    if (!addSeed(seedSet, seed3, 2, Merge()))
        addSeed(seedSet, seed3, Single());

    std::cout << "Merge Method:" << std::endl;
    for (TIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
        std::cout << "Seed: " << *it << std::endl;
    std::cout << std::endl;

    clear(seedSet);
    Score<int, Simple> scoreScheme(2, -1, -2);
    if (!addSeed(seedSet, seed0, 2, 1, scoreScheme, seqH, seqV, Chaos()))
        addSeed(seedSet, seed0, Single());
    if (!addSeed(seedSet, seed1, 2, 1, scoreScheme, seqH, seqV, Chaos()))
        addSeed(seedSet, seed1, Single());
    if (!addSeed(seedSet, seed2, 2, 1, scoreScheme, seqH, seqV, Chaos()))
        addSeed(seedSet, seed2, Single());
    if (!addSeed(seedSet, seed3, 2, 1, scoreScheme, seqH, seqV, Chaos()))
        addSeed(seedSet, seed3, Single());

    std::cout << "Chaos Method:" << std::endl;
    for (TIterator it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
        std::cout << "Seed: " << *it << std::endl;
    std::cout << std::endl;

}
