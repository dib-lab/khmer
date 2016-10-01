// A demonstration of using khmer to count k-mers; in this case, the request
// was for exact counting.

#include <vector>
#include <cmath>
#include "khmer.hh"
#include "counting.hh"

using namespace khmer;

int main()
{
    unsigned int ksize = 11;

    // For exact counting, you need to create one table that is >= 4**k. This
    // will be that size in bytes, note, so you will need the appropriate amount
    // of memory.

    // If `ksize` is even, note that k-mers will collapse with their reverse
    // complement. In that case, a table size of 4**(k-1) + k is required.

    std::vector<HashIntoType> tablesize;
    tablesize.push_back(pow(4, ksize));

    CountingHash ktable(ksize, tablesize);

    ktable.consume_string("ATGGCGATGGCAAGTAGGACCCAGATGGACCAAAG");

    std::cout << "count for: " << "ATGGCGATGGC" << " is " <<
        ktable.get_count("ATGGCGATGGC") << "\n";

    ktable.add("ATGGCGATGGC");

    std::cout << "count for: " << "ATGGCGATGGC" << " is " <<
        ktable.get_count("ATGGCGATGGC") << "\n";

    std::cout << "count for: " << "GTGGCGATGGC" << " is " <<
        ktable.get_count("GTGGCGATGGC") << "\n";

    return 0;
}
