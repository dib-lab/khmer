/***************************************************************************
 *  examples/containers/deque1.cpp
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
#include <stxxl/deque>
#include <iostream>

int main()
{
    typedef stxxl::deque<int> deque;
    deque my_deque;

    my_deque.push_front(2);
    my_deque.push_front(11);
    my_deque.push_back(5);
    my_deque.push_back(8);
    // deque now stores: |11|2|5|8|

    std::cout << "return 'first' element: " << my_deque.front() << std::endl; // prints 11
    std::cout << "return 'last' element: " << my_deque.back() << std::endl;   // prints 8
    std::cout << "random access: " << my_deque[2] << std::endl;               // prints 5

    // generate forward iterator
    stxxl::deque_iterator<deque> deque_iterator = my_deque.begin();

    // iterate over my_deque, access values and delete them afterwards
    while (!my_deque.empty())
    {
        std::cout << *deque_iterator << " ";
        ++deque_iterator;
        my_deque.pop_front();
    }

    return 0;
}
//! [example]
