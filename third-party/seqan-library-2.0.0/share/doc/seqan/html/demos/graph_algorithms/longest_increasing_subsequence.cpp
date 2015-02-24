#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    // Create a sequence of integers.
    String<unsigned> seq;
    appendValue(seq, 5);
    appendValue(seq, 3);
    appendValue(seq, 4);
    appendValue(seq, 9);
    appendValue(seq, 6);
    appendValue(seq, 2);
    appendValue(seq, 1);
    appendValue(seq, 8);
    appendValue(seq, 7);
    appendValue(seq, 10);

    // Compute the longest increasing subsequence.
    typedef Position<String<unsigned> >::Type TPosition;
    String<TPosition> pos;
    longestIncreasingSubsequence(seq, pos);

    // Print the result to stdout.
    for (unsigned i = 0; i < length(seq); ++i)
        std::cout << seq[i] << ',';

    std::cout << "\n"
              << "Lis: \n";
    for (int i = length(pos) - 1; i >= 0; --i)
        std::cout << seq[pos[i]] <<  ',';
    std::cout << "\n";

    return 0;
}
