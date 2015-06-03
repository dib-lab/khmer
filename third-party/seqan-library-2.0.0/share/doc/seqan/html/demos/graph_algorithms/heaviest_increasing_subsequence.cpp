#include <iostream>
#include <seqan/sequence.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main()
{
    // Fill a string and define corresponding weights.
    String<char> seq("zeitgeist");
    String<unsigned int> weights;
    resize(weights, length(seq), 1);
    assignProperty(weights, 2, 10);

    // Compute heaviest increasing subsequence.
    typedef Position<String<unsigned int> >::Type TPosition;
    String<TPosition> pos;
    unsigned int w = heaviestIncreasingSubsequence(seq, weights, pos);

    // Print the results to stdout.
    for (int i = 0; i < (int) length(seq); ++i)
        std::cout << seq[i] << "(Weight=" << getProperty(weights, i) << "),";
    std::cout << "\n"
              << "His: \n";
    for (int i = length(pos) - 1; i >= 0; --i)
        std::cout << seq[pos[i]] <<  ',';
    std::cout << "(Weight=" << w << ")\n";

    return 0;
}
