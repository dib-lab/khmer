#include <vector>
#include "khmer.hh"
#include "counting.hh"

using namespace khmer;

int main()
{
    std::vector<HashIntoType> tablesize;
    tablesize.push_back(4194304);
    
    CountingHash ktable(11, tablesize);

    ktable.consume_string("ATGGCGATGGCAAGTAGGACCCAGATGGACCAAAG");

    std::cout << "count for: " << "ATGGCGATGGC" << " is " <<
        ktable.get_count("ATGGCGATGGC") << "\n";

    ktable.count("ATGGCGATGGC");

    std::cout << "count for: " << "ATGGCGATGGC" << " is " <<
        ktable.get_count("ATGGCGATGGC") << "\n";

    std::cout << "count for: " << "GTGGCGATGGC" << " is " <<
        ktable.get_count("GTGGCGATGGC") << "\n";

    return 0;
}
