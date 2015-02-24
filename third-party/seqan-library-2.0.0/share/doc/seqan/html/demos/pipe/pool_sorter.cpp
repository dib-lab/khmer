#include <functional>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/pipe.h>

using namespace seqan;

// Custom 3-way comparator functor.  Return a value </==/> 0 (less than, equal
// to, greater than).
struct MyIntLessCmp :
    std::binary_function<int, int, int>
{
    int operator()(int const & lhs, int const & rhs) const
    {
        return lhs - rhs;
    }

};

int main()
{
    typedef MyIntLessCmp                     TLess;
    typedef SorterSpec<SorterConfig<TLess> > TSorterSpec;
    typedef Pool<int, TSorterSpec>           TSorterPool;

    // Fill the sorter pool.
    TSorterPool sorterPool;
    resize(sorterPool, 3);
    beginWrite(sorterPool);
    push(sorterPool, 3);
    push(sorterPool, 10);
    push(sorterPool, -1);
    endWrite(sorterPool);

    // Fetch the resorted elements from the pool and print to stdout.
    std::cout << "Sorted elements:\n";
    beginRead(sorterPool);
    for (; !eof(sorterPool); pop(sorterPool))
        std::cout << front(sorterPool) << "\n";
    endRead(sorterPool);

    return 0;
}
