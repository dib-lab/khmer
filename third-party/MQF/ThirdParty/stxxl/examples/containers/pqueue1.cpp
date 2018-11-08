/***************************************************************************
 *  examples/containers/pqueue1.cpp
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
#include <stxxl/priority_queue>
#include <iostream>
#include <limits>

// comparison struct for priority queue where top() returns the smallest contained value:
struct ComparatorGreater
{
    bool operator () (const int& a, const int& b) const
    { return (a > b); }

    int min_value() const
    { return std::numeric_limits<int>::max(); }
};

int main()
{
    typedef stxxl::PRIORITY_QUEUE_GENERATOR<int, ComparatorGreater, 128*1024*1024, 1024*1024>::result pqueue_type;
    typedef pqueue_type::block_type block_type;

    // block_type::raw_size = 262144 bytes
    // use 64 block read and write pools each to enable overlapping between I/O and computation
    const unsigned int mem_for_pools = 32 * 1024 * 1024;
    stxxl::read_write_pool<block_type> pool((mem_for_pools / 2) / block_type::raw_size, (mem_for_pools / 2) / block_type::raw_size);
    pqueue_type my_pqueue(pool);  // creates stxxl priority queue instance with read-write-pool

    my_pqueue.push(5);
    my_pqueue.push(4);
    my_pqueue.push(19);
    my_pqueue.push(1);
    assert(my_pqueue.size() == 4);

    assert(my_pqueue.top() == 1);
    STXXL_MSG("Smallest inserted value in my: " << my_pqueue.top());

    my_pqueue.pop();  // pop the 1 on top

    assert(my_pqueue.top() == 4);
    STXXL_MSG("Smallest value after 1 pop(): " << my_pqueue.top());

    return 0;
}
//! [example]
