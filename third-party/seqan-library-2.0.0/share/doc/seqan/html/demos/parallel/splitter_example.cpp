#include <iostream>
#include <seqan/parallel.h>

using namespace seqan;

int main()
{
    // instantiate a Splitter to divide the interval [10,20) in 3 subintervals
    Splitter<unsigned> splitter(10, 20, 3);

    // output all subintervals
    for (unsigned i = 0; i < length(splitter); ++i)
        std::cout << '[' << splitter[i] << ',' << splitter[i + 1] << ')' << std::endl;

    return 0;
}
