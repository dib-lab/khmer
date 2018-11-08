/***************************************************************************
 *  examples/containers/queue2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/queue>
#include <iostream>

int main()
{
    // template parameter <value_type, block_size, allocation_strategy, size_type>
    typedef stxxl::queue<unsigned int> a_queue;

    // construct queue with default parameters
    a_queue my_queue;

    unsigned int random;
    stxxl::random_number32 rand32;  // define random number generator
    stxxl::uint64 number_of_elements = 64 * 1024 * 1024;

    // push random values in the queue
    for (stxxl::uint64 i = 0; i < number_of_elements; i++)
    {
        random = rand32();  // generate random integers from interval [0,2^32)
        my_queue.push(random);
    }

    unsigned int last_inserted = my_queue.back();
    STXXL_MSG("last element inserted: " << last_inserted);

    // identify smaller element than first_inserted, search in growth-direction (front->back)
    while (!my_queue.empty())
    {
        if (last_inserted > my_queue.front())
        {
            STXXL_MSG("found smaller element: " << my_queue.front() << " than last inserted element");
            break;
        }
        std::cout << my_queue.front() << " " << std::endl;
        my_queue.pop();
    }

    return 0;
}
