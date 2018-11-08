/***************************************************************************
 *  tests/containers/test_queue2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#define STXXL_VERBOSE_LEVEL 0

// stxxl::queue contains deprecated funtions
#define STXXL_NO_DEPRECATED 1

#include <stxxl/queue>

typedef stxxl::uint64 my_type;

// forced instantiation
template class stxxl::queue<my_type>;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " [n in MiB]" << std::endl;
        return -1;
    }

    stxxl::int64 mega = 1 << 20;
    stxxl::int64 megabytes = atoi(argv[1]);
    stxxl::int64 nelements = megabytes * mega / sizeof(my_type);
    stxxl::int64 i;
    my_type in = 0, out = 0;

    stxxl::queue<my_type> q;

    STXXL_MSG("op-sequence: ( push, pop, push ) * n");
    for (i = 0; i < nelements; ++i)
    {
        if ((i % mega) == 0)
            STXXL_MSG("Insert " << i);

        q.push(in++);

        STXXL_CHECK(q.front() == out);
        q.pop();
        ++out;

        q.push(in++);
    }
    STXXL_MSG("op-sequence: ( pop, push, pop ) * n");
    for ( ; i > 0; --i)
    {
        if ((i % mega) == 0)
            STXXL_MSG("Remove " << i);

        STXXL_CHECK(q.front() == out);
        q.pop();
        ++out;

        q.push(in++);

        STXXL_CHECK(q.front() == out);
        q.pop();
        ++out;
    }
    STXXL_CHECK(q.empty());
    STXXL_CHECK(in == out);

    std::cout << *stxxl::stats::get_instance();

    return 0;
}
