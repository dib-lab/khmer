/***************************************************************************
 *  examples/containers/sorter1.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! [example]
#include <stxxl/sorter>
#include <stxxl/stats>
#include <stxxl/timer>
#include <limits>

struct my_comparator
{
    bool operator () (const int& a, const int& b) const
    {
        return a < b;
    }

    int min_value() const
    {
        return std::numeric_limits<int>::min();
    }

    int max_value() const
    {
        return std::numeric_limits<int>::max();
    }
};

int main()
{
    // template parameter <ValueType, CompareType, BlockSize, AllocStr(optional)>
    typedef stxxl::sorter<int, my_comparator> sorter_type;

    // create sorter object (CompareType(), MainMemoryLimit)
    sorter_type int_sorter(my_comparator(), 64 * 1024 * 1024);

    // fill sorter with elements order in descending order
    for (int i = 1000; i > 0; i--)
    {
        int_sorter.push(i);
    }

    int_sorter.sort();  // sort elements (in ascending order)

    // walk through sorted values and print them out
    while (!int_sorter.empty())
    {
        std::cout << *int_sorter << " ";
        ++int_sorter;
    }

    return 0;
}
//! [example]
