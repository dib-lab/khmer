/***************************************************************************
 *  tests/containers/test_sequence.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iterator>
#include <stxxl/sequence>
#include <stxxl/random>

int main(int argc, char* argv[])
{
    stxxl::uint64 ops = (argc >= 2) ? stxxl::atouint64(argv[1]) : 32 * 1024 * 1024;

    stxxl::random_number32 random;
    stxxl::sequence<int> XXLDeque;
    std::deque<int> STDDeque;

    for (stxxl::uint64 i = 0; i < ops; ++i)
    {
        unsigned curOP = random() % 6;
        unsigned value = random();
        switch (curOP)
        {
        case 0: // make insertion a bit more likely
        case 1:
            XXLDeque.push_front(value);
            STDDeque.push_front(value);
            break;
        case 2: // make insertion a bit more likely
        case 3:
            XXLDeque.push_back(value);
            STDDeque.push_back(value);
            break;
        case 4:
            if (!XXLDeque.empty())
            {
                XXLDeque.pop_front();
                STDDeque.pop_front();
            }
            break;
        case 5:
            if (!XXLDeque.empty())
            {
                XXLDeque.pop_back();
                STDDeque.pop_back();
            }
            break;
        }

        STXXL_CHECK(XXLDeque.empty() == STDDeque.empty());
        STXXL_CHECK(XXLDeque.size() == STDDeque.size());

        if (XXLDeque.size() > 0)
        {
            STXXL_CHECK(XXLDeque.back() == STDDeque.back());
            STXXL_CHECK(XXLDeque.front() == STDDeque.front());
        }

        if (!(i % 1000000))
        {
            std::cout << "Complete check of sequence/deque (size " << XXLDeque.size() << ")\n";
            stxxl::sequence<int>::stream stream = XXLDeque.get_stream();
            std::deque<int>::const_iterator b = STDDeque.begin();

            while (!stream.empty())
            {
                STXXL_CHECK(b != STDDeque.end());
                STXXL_CHECK(*stream == *b);
                ++stream;
                ++b;
            }

            STXXL_CHECK(b == STDDeque.end());
        }

        if (!(i % 1000000))
        {
            std::cout << "Complete check of reverse sequence/deque (size " << XXLDeque.size() << ")\n";
            stxxl::sequence<int>::reverse_stream stream = XXLDeque.get_reverse_stream();
            std::deque<int>::reverse_iterator b = STDDeque.rbegin();

            while (!stream.empty())
            {
                STXXL_CHECK(b != STDDeque.rend());
                STXXL_CHECK(*stream == *b);
                ++stream;
                ++b;
            }

            STXXL_CHECK(b == STDDeque.rend());
        }
    }

    return 0;
}
