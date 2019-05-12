/***************************************************************************
 *  examples/containers/deque2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/deque>
#include <iostream>

int main()
{
    typedef stxxl::deque<unsigned int> deque;
    deque my_deque;

    unsigned int random, p, x;
    unsigned int smaller_left = 0;
    unsigned int smaller_right = 0;
    stxxl::random_number32 rand32;
    stxxl::uint64 number_of_elements = 8 * 1024 * 1024;

    // fill deque with random integer values
    for (stxxl::uint64 i = 0; i < number_of_elements; i++)
    {
        random = rand32();  // produce random integer from intervall [0,2^32)
        my_deque.push_front(random);
    }

    stxxl::deque_iterator<deque> deque_iterator = my_deque.begin();

    // Access random element x at position p(x) in the deque
    p = (unsigned int)(rand32() % number_of_elements);
    x = my_deque[p];

    // Count number of smaller elements from the front to p(x) - 1
    for (stxxl::uint64 j = 0; j < p; j++)
    {
        if (*deque_iterator < x)
        {
            smaller_left += 1;
        }
        ++deque_iterator;
    }

    ++deque_iterator;

    // Count number of smaller elements from p(x) + 1 to the end
    for (stxxl::uint64 k = p + 1; k < number_of_elements - 1; k++)
    {
        if (*deque_iterator < x)
        {
            smaller_right += 1;
        }
        ++deque_iterator;
    }

    STXXL_MSG("smaller left: " << smaller_left << ", smaller right: " << smaller_right);

    return 0;
}
