// A demonstration of using khmer for exact k-mer counting. The memory required
// is 4^k, which limits this to small values of k.

#include <vector>
#include <cmath>
#include "oxli/oxli.hh"
#include "oxli/hashgraph.hh"


using namespace oxli;

int main()
{
    unsigned int ksize = 11;
    uint64_t nkmers = pow(4, ksize);

    // For exact counting, you need to create one table that is >= 4**k. This
    // will be that size in bytes, note, so you will need the appropriate amount
    // of memory.

    // If `ksize` is even, note that k-mers will collapse with their reverse
    // complement. In that case, a table size of 4**(k-1) + k is required.

    std::vector<uint64_t> tablesize = {nkmers};
    Countgraph counts(ksize, tablesize);

    counts.consume_string("ATGGCGATGGCAAGTAGGACCCAGATGGACCAAAG");

    std::cout << "count for: " << "ATGGCGATGGC" << " is " <<
        counts.get_count("ATGGCGATGGC") << "\n";

    counts.add("ATGGCGATGGC");

    std::cout << "count for: " << "ATGGCGATGGC" << " is " <<
        counts.get_count("ATGGCGATGGC") << "\n";

    std::cout << "count for: " << "GTGGCGATGGC" << " is " <<
        counts.get_count("GTGGCGATGGC") << "\n";

    return 0;
}
