// A demonstration of khmer's primary sequence loading function.

#include <vector>
#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"

using namespace oxli;
using namespace read_parsers;

int main(int argc, char **argv)
{
    unsigned int ksize = 21;

    // Initialize a Count-min sketch with 4 hash functions (4 distinct tables
    // with a prime number of buckets). The sum of these values is the memory
    // consumption of the Bloom filter in bits. See `khmer.get_n_primes_near_x`
    // from the Python API.
    std::vector<uint64_t> tablesizes = {
        499999897, 499999909, 499999931, 499999993
    };
    Counttable counts(ksize, tablesizes);

    // Read the input file into the Count-min sketch
    unsigned int seqs_consumed = 0;
    unsigned long long kmers_consumed = 0;
    counts.consume_seqfile<FastxReader>(argv[1], seqs_consumed, kmers_consumed);
    std::cout << "Loaded " << seqs_consumed << " sequences and "
              << kmers_consumed << " k-mers from " << argv[1] << '\n';

    // Do some k-mer abundance queries
    std::cout << "The k-mer 'CAGCGCCGTGTTGTTGCAATT' appears "
              << counts.get_count("CAGCGCCGTGTTGTTGCAATT")
              << " times in the data.\n";
    std::cout << "The k-mer 'GATTACAGATTACAGATTACA' appears "
              << counts.get_count("GATTACAGATTACAGATTACA")
              << " times in the data.\n";

    return 0;
}
