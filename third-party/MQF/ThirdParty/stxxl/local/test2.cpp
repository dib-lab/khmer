/***************************************************************************
 *  local/test2.cpp
 *
 *  This is another example file included in the local/ directory of STXXL. All
 *  .cpp files in local/ are automatically compiled and linked with STXXL by
 *  CMake.  You can use this method for simple prototype applications.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <limits>

#include <stxxl/random>
#include <stxxl/sorter>

struct my_less : std::less<int64_t>
{
    int64_t min_value() const { return std::numeric_limits<int64_t>::min(); }
    int64_t max_value() const { return std::numeric_limits<int64_t>::max(); }
};

int main()
{
    stxxl::scoped_print_timer
        timer("overall work", 600 * 1024 * 1024 * (int64_t)sizeof(int64_t));

    // create sorter
    stxxl::sorter<int64_t, my_less> sorter(my_less(), 256 * 1024 * 1024);

    // fill sorter with random integers
    {
        stxxl::scoped_print_timer
            timer("presort+write random numbers", 600 * 1024 * 1024 * (int64_t)sizeof(int64_t));

        stxxl::random_number32 random;

        for (size_t i = 0; i < 600 * 1024 * 1024; ++i) {
            sorter.push(random() * random());
        }
    }

    sorter.sort();

    // get data back in sorted order
    {
        stxxl::scoped_print_timer
            timer("read+merge random numbers", 600 * 1024 * 1024 * (int64_t)sizeof(int64_t));

        int64_t first = *sorter, last = first, count = 1;
        ++sorter;

        while (!sorter.empty())
            last = *sorter, ++sorter, ++count;

        // output first and last items:
        std::cout << count << " items sorted ranging from "
                  << first << " to " << last << std::endl;
    }

    return 0;
}
