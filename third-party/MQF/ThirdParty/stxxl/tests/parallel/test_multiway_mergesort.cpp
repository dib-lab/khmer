/***************************************************************************
 *  tests/parallel/test_multiway_mergesort.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/parallel.h>
#include <stxxl/bits/parallel/multiway_mergesort.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/is_sorted.h>
#include <stxxl/random>
#include <iostream>

struct Something
{
    int a, b;

    Something(int x = 0)
        : a(x), b(x * x)
    { }

    bool operator < (const Something& other) const
    {
        return a < other.a;
    }

    friend std::ostream& operator << (std::ostream& os, const Something& s)
    {
        return os << '(' << s.a << ',' << s.b << ')';
    }
};

template <bool Stable>
void test_size(unsigned int size)
{
    std::cout << "testing parallel_sort_mwms with " << size << " items.\n";

    std::vector<Something> v(size);
    std::less<Something> cmp;

    stxxl::random_number32 rnd;

    for (unsigned int i = 0; i < size; ++i)
        v[i] = rnd();

    stxxl::parallel::parallel_sort_mwms<Stable>(v.begin(), v.end(), cmp, 8);

    STXXL_CHECK(stxxl::is_sorted(v.begin(), v.end(), cmp));
}

int main()
{
    // run multiway mergesort tests for 0..256 sequences
    for (unsigned int i = 0; i < 256; ++i)
    {
        test_size<false>(i);
        test_size<true>(i);
    }

    // run multiway mergesort tests for 0..256 sequences
    for (unsigned int i = 256; i <= 16 * 1024 * 1024; i = 2 * i - i / 2)
    {
        test_size<false>(i);
        test_size<true>(i);
    }

    return 0;
}
