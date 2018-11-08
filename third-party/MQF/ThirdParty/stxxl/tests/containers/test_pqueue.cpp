/***************************************************************************
 *  tests/containers/test_pqueue.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/test_pqueue.cpp
//! This is an example of how to use \c stxxl::PRIORITY_QUEUE_GENERATOR
//! and \c stxxl::priority_queue

#include <limits>
#include <stxxl/priority_queue>
#include <stxxl/timer>

using stxxl::uint64;
using stxxl::scoped_print_timer;

#define RECORD_SIZE 128

const uint64 volume = 1024 * 1024; // in KiB

struct my_type
{
    typedef int key_type;
    key_type key;
    char data[RECORD_SIZE - sizeof(key_type)];
    my_type() { }
    explicit my_type(key_type k) : key(k)
    {
#if STXXL_WITH_VALGRIND
        memset(data, 0, sizeof(data));
#endif
    }

    friend std::ostream& operator << (std::ostream& o, const my_type& obj)
    {
        return o << obj.key;
    }
};

struct my_cmp : std::binary_function<my_type, my_type, bool> // greater
{
    bool operator () (const my_type& a, const my_type& b) const
    {
        return a.key > b.key;
    }

    my_type min_value() const
    {
        return my_type(std::numeric_limits<my_type::key_type>::max());
    }
};

// forced instantiation
template class stxxl::PRIORITY_QUEUE_GENERATOR<
        my_type, my_cmp, 32* 1024* 1024, volume / sizeof(my_type)
        >;

int main()
{
/*
  unsigned BufferSize1_ = 32, // equalize procedure call overheads etc.
  unsigned N_ = 512, // bandwidth
  unsigned IntKMAX_ = 64, // maximal arity for internal mergers
  unsigned IntLevels_ = 4,
  unsigned BlockSize_ = (2*1024*1024),
  unsigned ExtKMAX_ = 64, // maximal arity for external mergers
  unsigned ExtLevels_ = 2,
 */
    //typedef priority_queue<priority_queue_config<my_type,my_cmp,
    //  32,512,64,3,(4*1024),0x7fffffff,1> > pq_type;

    typedef stxxl::PRIORITY_QUEUE_GENERATOR<
            my_type, my_cmp, 32* 1024* 1024, volume / sizeof(my_type)
            > gen;
    typedef gen::result pq_type;
    typedef pq_type::block_type block_type;

    STXXL_MSG("Block size: " << block_type::raw_size);
    STXXL_MSG("AI: " << gen::AI);
    STXXL_MSG("X : " << gen::X);
    STXXL_MSG("N : " << gen::N);
    STXXL_MSG("AE: " << gen::AE);

    const unsigned mem_for_pools = 128 * 1024 * 1024;
    stxxl::read_write_pool<block_type> pool(
        (mem_for_pools / 2) / block_type::raw_size,
        (mem_for_pools / 2) / block_type::raw_size
        );
    pq_type p(pool);

    uint64 nelements = volume * 1024 / sizeof(my_type);
    STXXL_MSG("Internal memory consumption of the priority queue: " << p.mem_cons() << " B");
    STXXL_MSG("Max elements: " << nelements);

    {
        scoped_print_timer timer("Filling PQ",
                                 nelements * sizeof(my_type));

        for (uint64 i = 0; i < nelements; i++)
        {
            if ((i % (1024 * 1024)) == 0)
                STXXL_MSG("Inserting element " << i);
            p.push(my_type(int(nelements - i)));
        }
    }

    STXXL_CHECK(p.size() == (stxxl::uint64)nelements);

#if 0
    // test swap
    pq_type p1(pool);
    std::swap(p, p1);
    std::swap(p, p1);
#endif

    STXXL_MSG("Internal memory consumption of the priority queue: " << p.mem_cons() << " B");

    {
        scoped_print_timer timer("Emptying PPQ",
                                 nelements * sizeof(my_type));

        for (uint64 i = 0; i < nelements; ++i)
        {
            STXXL_CHECK(!p.empty());
            //STXXL_MSG( p.top() );
            STXXL_CHECK(p.top().key == int(i + 1));
            p.pop();
            if ((i % (1024 * 1024)) == 0)
                STXXL_MSG("Element " << i << " popped");
        }
    }

    STXXL_CHECK(p.size() == 0);
    STXXL_CHECK(p.empty());

    STXXL_MSG("Internal memory consumption of the priority queue: " << p.mem_cons() << " B");

    return 0;
}
