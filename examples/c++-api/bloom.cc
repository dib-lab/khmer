// A demonstration of using khmer to query a dataset for a k-mer. Typically
// khmer accrues a small false positive rate in order to save substantially on
// memory requirements.

#include <vector>
#include "khmer.hh"
#include "hashtable.hh"

using namespace khmer;

int main()
{
    unsigned int ksize = 21;

    // Initialize a Bloom filter with 4 hash functions (4 distinct tables with a
    // prime number of buckets). The sum of these values is the memory
    // consumption of the Bloom filter in bits. See `khmer.get_n_primes_near_x`
    // from the Python API.
    std::vector<uint64_t> tablesizes = {
        499999897, 499999909, 499999931, 499999993
    };
    Nodetable ktable(ksize, tablesizes);

    ktable.consume_string("GCTGCACCGATGTACGCAAAGCTATTTAAAACCATAACTATTCTCACTTA");

    std::cout << "count for: " << "GCTGCACCGATGTACGCAAAG" << " is " <<
        ktable.get_count("GCTGCACCGATGTACGCAAAG") << "\n";

    ktable.add("GCTGCACCGATGTACGCAAAG");

    std::cout << "count for: " << "GCTGCACCGATGTACGCAAAG" << " is " <<
        ktable.get_count("GCTGCACCGATGTACGCAAAG") << "\n";

    std::cout << "count for: " << "GATTACAGATTACAGATTACA" << " is " <<
        ktable.get_count("GATTACAGATTACAGATTACA") << "\n";

    return 0;
}
