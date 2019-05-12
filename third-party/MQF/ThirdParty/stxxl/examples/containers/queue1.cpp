/***************************************************************************
 *  examples/containers/queue1.cpp
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
#include <stxxl/queue>
#include <iostream>

int main()
{
    typedef stxxl::queue<unsigned int> queue;
    queue my_queue;

    my_queue.push(5);
    my_queue.push(11);
    my_queue.push(3);
    my_queue.push(7);
    // my_queue now stores: |7|3|11|5|

    assert(my_queue.size() == 4);

    std::cout << "back element " << my_queue.back() << std::endl;   // prints out 7 (last inserted element)
    assert(my_queue.back() == 7);
    std::cout << "front element " << my_queue.front() << std::endl; // prints out 5 (first inserted element)
    assert(my_queue.front() == 5);
    my_queue.pop();                                                 // deletes element 5, queue now stores: |7|3|11|

    std::cout << "front element " << my_queue.front() << std::endl; // prints out 11
    assert(my_queue.front() == 11);

    return 0;
}
//! [example]
