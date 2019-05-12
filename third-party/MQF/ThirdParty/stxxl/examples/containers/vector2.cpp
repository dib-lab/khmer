/***************************************************************************
 *  examples/containers/vector2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/vector>
#include <iostream>

int main()
{
    // template parameter <value_type, page_size, number_of_pages, block_size, alloc_strategy, paging_strategy>
    typedef stxxl::VECTOR_GENERATOR<unsigned int, 4, 8, 1*1024*1024, stxxl::RC, stxxl::lru>::result vector_type;

    vector_type my_vector;
    unsigned int counter = 0;
    unsigned int tmp;
    stxxl::random_number<> rand;
    stxxl::uint64 number_of_elements = 32 * 1024 * 1024;

    // fill vector with random integers
    for (stxxl::uint64 i = 0; i < number_of_elements; ++i)
    {
        tmp = rand(123456789);  // generate random number from the interval [0,123456789)
        my_vector.push_back(tmp);
    }

    // construct iterator
    vector_type::const_iterator iter = my_vector.begin();

    // use iterator to advance my_vector and calculate number of even elements
    for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
    {
        //std::cout << *iter << " ";
        if (*iter % 2 == 0)  // is my_vector's current element even?
        {
            ++counter;
        }
        iter++;
    }

    STXXL_MSG("found " << counter << " even numbers in V");
    return 0;
}
