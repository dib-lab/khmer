///A tutorial about suffix arrays.
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<char> text = "hello world!";
    String<char> pattern = "l";
    String<unsigned> sa;

//Build a suffix array using the Skew7 algorithm.
    resize(sa, length(text));
    createSuffixArray(sa, text, Skew7());

//Search the interval of suffices beginning with the pattern.
    Pair<unsigned> hitRange;
    hitRange = equalRangeSA(text, sa, pattern);

//Output the suffix indices, i.e. the occurrences of the pattern.
    for (unsigned i = hitRange.i1; i < hitRange.i2; ++i)
        std::cout << sa[i] << " ";
    std::cout << std::endl;

    return 0;
}
