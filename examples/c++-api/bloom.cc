// A demonstration of using khmer to query a dataset for a k-mer. Typically
// khmer accrues a small false positive rate in order to save substantially on
// memory requirements.

#include <vector>
#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"

using namespace oxli;

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
    Nodetable bloomfilter(ksize, tablesizes);

    bloomfilter.consume_string("GCTGCACCGATGTACGCAAAGCTATTTAAAACCATAACTATTCTCACTTA");

    std::cout << "count for: 'GCTGCACCGATGTACGCAAAG' is "
              << bloomfilter.get_count("GCTGCACCGATGTACGCAAAG") << "\n";

    bloomfilter.add("GCTGCACCGATGTACGCAAAG");

    std::cout << "count for: 'GCTGCACCGATGTACGCAAAG' is "
              << bloomfilter.get_count("GCTGCACCGATGTACGCAAAG") << "\n";

    std::cout << "count for: 'GATTACAGATTACAGATTACA' is "
              << bloomfilter.get_count("GATTACAGATTACAGATTACA") << "\n";

    return 0;
}
