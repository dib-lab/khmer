#include <string>
#include <memory>
#include <cassert>
#include <iostream>

#include "cyclichash.h"

/**
* Test of the prepend and append functions to test slightly longer and slightly shorter n-grams.
*/

int main(int argc, char * argv[])
{
    CyclicHash<uint64_t>  hf(4, 64);
    string input = "XABCDY";
    string base(input.begin() + 1, input.end() - 1);
    string extend(input.begin() + 1, input.end());
    string prepend(input.begin(), input.end() - 1);

    for (string::const_iterator j = base.begin(); j != base.end(); ++j)
    {
        hf.eat(*j);
    }

    std::cout << base << " "  << hf.hash(base) << std::endl;
    std::cout << prepend << " " << hf.hash_prepend(input[0]) << " " << hf.hash(prepend) << std::endl;
    std::cout << extend << " " << hf.hash_extend(input.back()) << " " << hf.hash(extend) << std::endl;

    assert(hf.hashvalue == hf.hash(base));
    assert(hf.hash_prepend(input[0]) == hf.hash(prepend));
    assert(hf.hash_extend(input.back()) == hf.hash(extend));

    return 0;

}
