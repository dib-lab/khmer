/***************************************************************************
 *  tests/containers/test_queue.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2005 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

// stxxl::queue contains deprecated funtions
#define STXXL_NO_DEPRECATED 1

#include <queue>
#include <stxxl/queue>

typedef unsigned my_type;

template <class q1type, class q2type>
void check(const q1type& q1, const q2type& q2)
{
    STXXL_CHECK(q1.empty() == q2.empty());
    STXXL_CHECK(q1.size() == q2.size());
    if (!q1.empty())
    {
        if (q1.front() != q2.front() || q1.back() != q2.back())
            STXXL_MSG(q1.size() << ": (" << q1.front() << ", " << q1.back() << ") (" << q2.front() << ", " << q2.back() << ")" << (q1.front() == q2.front() ? "" : " FRONT"));
        STXXL_CHECK(q1.front() == q2.front());
        STXXL_CHECK(q1.back() == q2.back());
    }
}

// forced instantiation
template class stxxl::queue<my_type>;

int main()
{
#if 1
    //stxxl::set_seed(424648671);  // works fine with a single disk, fails with two disks
    STXXL_MSG("SEED=" << stxxl::get_next_seed());
    stxxl::srandom_number32(stxxl::get_next_seed());
    stxxl::random_number32_r rnd;
#else
    //stxxl::ran32State = 1028675152; // fails with two disks
    STXXL_MSG("ran32State=" << stxxl::ran32State);
    stxxl::random_number32 rnd;
#endif

    unsigned cnt;
    STXXL_MSG("Elements in a block: " << stxxl::queue<my_type>::block_type::size);

    // FIXME: can this be raised to recommended (3, 2) without breaking "Testing special case 4" or any other tests?
    stxxl::queue<my_type> xqueue(2, 2, -1);
    std::queue<my_type> squeue;
    check(xqueue, squeue);

    STXXL_MSG("Testing special case 4");
    cnt = stxxl::queue<my_type>::block_type::size;
    while (cnt--)
    {
        my_type val = rnd();
        xqueue.push(val);
        squeue.push(val);
        check(xqueue, squeue);
    }
    cnt = stxxl::queue<my_type>::block_type::size;
    while (cnt--)
    {
        xqueue.pop();
        squeue.pop();
        check(xqueue, squeue);
    }
    STXXL_MSG("Testing other cases ");
    cnt = 100 * stxxl::queue<my_type>::block_type::size;
    while (cnt--)
    {
        if ((cnt % 1000000) == 0)
            STXXL_MSG("Operations left: " << cnt << " queue size " << squeue.size());

        int rndtmp = rnd() % 3;
        if (rndtmp > 0 || squeue.empty())
        {
            my_type val = rnd();
            xqueue.push(val);
            squeue.push(val);
            check(xqueue, squeue);
        }
        else
        {
            xqueue.pop();
            squeue.pop();
            check(xqueue, squeue);
        }
    }
    while (!squeue.empty())
    {
        if ((cnt++ % 1000000) == 0)
            STXXL_MSG("Operations: " << cnt << " queue size " << squeue.size());

        int rndtmp = rnd() % 4;
        if (rndtmp >= 3)
        {
            my_type val = rnd();
            xqueue.push(val);
            squeue.push(val);
            check(xqueue, squeue);
        }
        else
        {
            xqueue.pop();
            squeue.pop();
            check(xqueue, squeue);
        }
    }
    cnt = 10 * stxxl::queue<my_type>::block_type::size;
    while (cnt--)
    {
        if ((cnt % 1000000) == 0)
            STXXL_MSG("Operations left: " << cnt << " queue size " << squeue.size());

        int rndtmp = rnd() % 3;
        if (rndtmp > 0 || squeue.empty())
        {
            my_type val = rnd();
            xqueue.push(val);
            squeue.push(val);
            check(xqueue, squeue);
        }
        else
        {
            xqueue.pop();
            squeue.pop();
            check(xqueue, squeue);
        }
    }
    typedef stxxl::queue<my_type>::block_type block_type;
    stxxl::read_write_pool<block_type> pool(5, 5);
    stxxl::queue<my_type> xqueue1(pool, -1);
    std::queue<my_type> squeue1;

    cnt = 10 * stxxl::queue<my_type>::block_type::size;
    while (cnt--)
    {
        if ((cnt % 1000000) == 0)
            STXXL_MSG("Operations left: " << cnt << " queue size " << squeue1.size());

        int rndtmp = rnd() % 3;
        if (rndtmp > 0 || squeue1.empty())
        {
            my_type val = rnd();
            xqueue1.push(val);
            squeue1.push(val);
            check(xqueue1, squeue1);
        }
        else
        {
            xqueue1.pop();
            squeue1.pop();
            check(xqueue1, squeue1);
        }
    }

    {
        // test proper destruction of a single-block queue
        stxxl::queue<int> q;
        q.push(42);
    }

    return 0;
}
