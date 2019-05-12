/***************************************************************************
 *  examples/containers/pqueue2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/priority_queue>
#include <limits>

// comparison struct for priority queue where top() returns the biggest contained value:
struct Cmp
{
    bool operator () (const int& a, const int& b) const
    { return a < b; }

    int min_value() const
    { return std::numeric_limits<int>::min(); }
};

int main()
{
    // use 64 GiB on main memory and 1 billion items at most
    typedef stxxl::PRIORITY_QUEUE_GENERATOR<int, Cmp, 64*1024*1024, 1024*1024>::result pq_type;
    typedef pq_type::block_type block_type;

    // block_type::raw_size = 262144 bytes
    // use 64 block read and write pools each to enable overlapping between I/O and computation
    const unsigned int mem_for_pools = 32 * 1024 * 1024;
    stxxl::read_write_pool<block_type> pool((mem_for_pools / 2) / block_type::raw_size, (mem_for_pools / 2) / block_type::raw_size);

    pq_type Q(pool);

    Q.push(1);
    Q.push(4);
    Q.push(2);
    Q.push(8);
    Q.push(5);
    Q.push(7);

    assert(Q.size() == 6);

    assert(Q.top() == 8);
    Q.pop();

    assert(Q.top() == 7);
    Q.pop();

    assert(Q.top() == 5);
    Q.pop();

    assert(Q.top() == 4);
    Q.pop();

    assert(Q.top() == 2);
    Q.pop();

    assert(Q.top() == 1);
    Q.pop();

    assert(Q.empty());

    return 0;
}
