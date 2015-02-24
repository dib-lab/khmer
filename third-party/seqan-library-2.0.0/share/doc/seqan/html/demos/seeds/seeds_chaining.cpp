#include <iostream>
#include <seqan/seeds.h>

using namespace seqan;

int main()
{
    // Build SeedSet.
    SeedSet<Seed<Simple>, Unordered> seedSet;
    addSeed(seedSet, Seed<Simple>(0, 93, 281, 342), Single());
    addSeed(seedSet, Seed<Simple>(3, 237, 127, 364), Single());
    addSeed(seedSet, Seed<Simple>(3, 284, 86, 368), Single());
    addSeed(seedSet, Seed<Simple>(5, 146, 239, 374), Single());
    addSeed(seedSet, Seed<Simple>(299, 352, 405, 460), Single());

    // Perform sparse chaining, uses time O(n log n).
    String<Seed<Simple> > chain;
    chainSeedsGlobally(chain, seedSet, SparseChaining());

    // Print results to stdout.
    for (unsigned i = 0; i < length(chain); ++i)
        std::cout << "Seed(" << beginPositionH(chain[i]) << ", "
                  << beginPositionV(chain[i]) << ", " << endPositionH(chain[i])
                  << ", " << endPositionV(chain[i]) << ")\n";

    return 0;
}
